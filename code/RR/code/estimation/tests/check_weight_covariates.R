## 在 anchored_v1 运行：Rscript estimation/tests/check_weight_covariates.R
source("load_anchored.R", encoding="UTF-8")
source_anchored()

## 均匀的 64 个二元协变量格：选定变量的密度比权重有独立解析解。
dat <- expand.grid(id=1:10, x=0:1, z1=0:1,z2=0:1,z3=0:1,z4=0:1,z5=0:1,z6=0:1)
dat$y <- as.integer(dat$id <= 2 + dat$x + dat$z1 + dat$z4)
dat <- cbind(dat, rr_gop_scores(covariate_matrix(dat)))
ipd <- list(data=dat)
target <- setNames(c(.7,.3,.7,.3,.7,.3),paste0("z",1:6))
Z <- covariate_matrix(dat)
all <- get.weight(Z,target)
all.explicit <- get.weight(Z,unname(target))
stopifnot(identical(all,all.explicit))

for (cols in list(c("z3","z1"), "z2", character(0), names(target))) {
  wf <- get.weight(Z[,cols,drop=FALSE],target[cols])
  analytic <- rep(1,nrow(dat))
  for (col in cols) analytic <- analytic * ifelse(dat[[col]]==1,2*target[col],2*(1-target[col]))
  stopifnot(max(abs(wf$weight-analytic)) < 5e-5,
            abs(mean(wf$weight)-1) < 1e-12,
            identical(names(wf$ga), cols), identical(wf$covariates,cols))
  balance <- colSums(covariate_matrix(dat)*wf$weight)/sum(wf$weight)
  stopifnot(max(abs(balance[cols]-target[cols]),0) < 1e-5,
            max(abs(balance[setdiff(names(target),cols)]-.5),0) < 1e-5)
}
selected.Z <- Z[,c("z3","z1"),drop=FALSE]
wf <- get.weight(selected.Z,target[c("z3","z1")])
for (mu in list(target[c(6,3,2,1,5,4)], target[c("z3","z1")],
                unname(target[c("z3","z1")]))) {
  stopifnot(max(abs(get.weight(selected.Z,mu)$weight-wf$weight)) < 1e-12)
}
unnamed <- get.weight(unname(selected.Z),unname(target[c("z3","z1")]))
stopifnot(max(abs(unnamed$weight-wf$weight)) < 1e-12,
          is.null(names(unnamed$ga)), is.null(unnamed$covariates))
bad <- try(get.weight(Z[,"z2",drop=FALSE],c(z1=.7)),silent=TRUE)
bad.length <- try(get.weight(selected.Z,.7),silent=TRUE)
stopifnot(inherits(bad,"try-error"), inherits(bad.length,"try-error"))
## 列名不限定 z 前缀；使用有独立解析解的两列矩阵检查名称对齐。
other.Z <- cbind(first=dat$z1,second=dat$z2)
other <- get.weight(other.Z,c(second=.3,first=.7))
analytic <- ifelse(dat$z1==1,1.4,.6)*ifelse(dat$z2==1,.6,1.4)
stopifnot(max(abs(other$weight-analytic)) < 5e-5,
          identical(names(other$ga),colnames(other.Z)),
          identical(other$target_mean,c(first=.7,second=.3)))
cat("PASS: matrix inputs, analytic weights, named/unnamed alignment, single and zero columns\n")

## 非对角目标协方差，检查联合方差使用实际所选列，同时保留全维度结局模型。
target.cov <- diag(6)*.15 + matrix(.01,6,6)
dimnames(target.cov) <- list(names(target),names(target))
for (cols in list(c("z3","z1"), "z2", character(0))) {
  wf <- get.weight(Z[,cols,drop=FALSE],target[cols])
  for (param in c("RR","RD")) {
    fit <- suppressWarnings(get.estimate(op_model="linear", "IPD",ipd,wf$weight,wf$ga,
      target_mean=target,target_cov=target.cov,m_AD=500,param=param))
    stopifnot(ncol(fit$vb)==7L,length(fit$point.est)==8L,all(is.finite(fit$cov)))
    vfun <- if(param=="RR") var.mle.rr.joint else var.mle.rd.joint
    manual <- vfun(dat$y,dat$x,fit$point.est[1],fit$point.est[2:8],wf$ga,
      fit$va,fit$vb,as.matrix(dat[,cols,drop=FALSE]),target[cols],
      target_cov=target.cov[cols,cols,drop=FALSE],m_AD=500)
    idx <- length(cols)+seq_len(8)
    stopifnot(max(abs(fit$cov-manual[idx,idx,drop=FALSE])) < 1e-10)
  }
}
wf <- get.weight(Z[,"z2",drop=FALSE],target["z2"])
firth <- suppressWarnings(get.estimate(op_model="linear", "IPD",ipd,wf$weight,wf$ga,
  target_mean=target,target_cov=target.cov,m_AD=500,method="firth"))
stopifnot(all(is.finite(firth$cov)),ncol(firth$vb)==7L)
cat("PASS: RR/RD and Firth use selected weight design; outcome covariates remain intact\n")

## RLP/CMH 的全维度汇总输入与手动截取的目标汇总应给出同一方差。
cols <- c("z3","z1")
wf <- get.weight(Z[,cols,drop=FALSE],target[cols])
for (fun in list(fit_robust_lp_joint,get.weighted.RR.sandwich)) {
  full <- fun(dat,wf$weight,gamma=wf$ga,target_mean=target,target_cov=target.cov,
              m_AD=500,account_for_target_mean=TRUE)
  selected <- fun(dat,wf$weight,gamma=unname(wf$ga),target_mean=target[cols],
                  target_cov=target.cov[cols,cols],m_AD=500,
                  account_for_target_mean=TRUE,weight_covariates=cols)
  stopifnot(isTRUE(all.equal(full,selected,tolerance=1e-12)))
}
score.AD <- get.aggre(ipd,covariates="h_E")
score.wf <- get.weight(as.matrix(dat[,"h_E",drop=FALSE]),c(h_E=0))
stopifnot(length(score.AD$anchored.mean.x)==1L,
          abs(weighted.mean(dat$h_E,score.wf$weight)) < 1e-5)
cat("PASS: RLP/CMH covariance selection and optional derived-score calibration\n")
