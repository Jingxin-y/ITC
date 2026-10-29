## 在 anchored_v1 运行：Rscript data/tests/check_getProbRR_GOP.R
source("load_anchored.R", encoding="UTF-8")
source_anchored_data()

## 从任意已知三臂风险构造输入，再反解回来；不依赖风险表的实现。
p <- as.matrix(expand.grid(p_C=c(.02,.2,.55,.85),p_A=c(.03,.3,.7),p_B=c(.05,.4,.8)))
p <- rbind(p, c(p_C=1e-6,p_A=3e-6,p_B=7e-7))
logrr_AC <- log(p[,"p_A"]/p[,"p_C"])
logrr_BC <- log(p[,"p_B"]/p[,"p_C"])
loggop <- rowSums(qlogis(p))
recovered <- getProbRR_GOP(logrr_AC,logrr_BC,loggop)
stopifnot(identical(colnames(recovered),c("p_C","p_A","p_B")),
          max(abs(log(recovered/p))) < 1e-10)

## 两个 RR=1 时，三臂均为 expit(logGOP/3)，而不是两臂的 expit(logOP/2)。
g <- c(-12,-3,0,3,12)
equal <- getProbRR_GOP(0,0,g)
stopifnot(max(abs(equal-plogis(g/3))) < 1e-10,
          identical(dim(getProbRR_GOP(0,0,0)),c(1L,3L)),
          max(abs(getProbScalarRR_GOP(0,0,0)-.5)) < 1e-12)
swapped <- getProbRR_GOP(logrr_BC,logrr_AC,loggop)
stopifnot(max(abs(swapped-recovered[,c("p_C","p_B","p_A")])) < 1e-10)
cat("PASS: three-arm risk round trip, scalar/vector calls and arm symmetry\n")

## 暂时禁用两臂转换，证明新 DGM 不调用旧 getProbRR()/getProbScalarRR()。
old.vector <- get0("getProbRR", mode="function")
old.scalar <- get0("getProbScalarRR", mode="function")
getProbRR <- getProbScalarRR <- function(...) stop("Two-arm risk converter used by joint DGM")
make_design <- function(r) joint_rr_gop_design(.4,.2,.3,.5,-3,.2,.4,
                                               ess_ratio=r,direction=c(-1,1,-1,1,-1,1))
design1 <- make_design(.5)
design2 <- make_design(.8)
set.seed(12)
ad1 <- data.generation("AD",500,design1)
set.seed(12)
ad2 <- data.generation("AD",500,design2)
set.seed(12)
ipd1 <- data.generation("IPD",500,design1)
stopifnot(identical(design1$e,design2$e),!identical(design1$q,design2$q),
          identical(design1$risk_table,design2$risk_table),
          identical(ad1$data,ad2$data),identical(ad1$truth,ad2$truth),
          identical(ad1$design$covariate_prob,design1$e),
          identical(ipd1$design$covariate_prob,design1$q),
          all((design1$q-design1$e)*c(-1,1,-1,1,-1,1)>0))
getProbRR <- old.vector
getProbScalarRR <- old.scalar
cat("PASS: AD stays fixed, IPD varies with ESS/direction, joint DGM excludes two-arm getProb\n")
