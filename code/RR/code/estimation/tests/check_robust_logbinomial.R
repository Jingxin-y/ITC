## 在 anchored_v1 根目录运行此文件。
source("load_anchored.R", encoding="UTF-8")
source_anchored()

dat <- expand.grid(id=1:30,x=0:1,z1=0:1,z2=0:1,z3=0:1)
dat$y <- as.integer(dat$id <= 3+2*dat$x+4*dat$z1+2*dat$z2+3*dat$z3+
                     2*dat$x*dat$z1+dat$x*dat$z3)
cols <- c("z3","z1")
target <- c(z1=.4,z2=.5,z3=.6)
Sigma <- matrix(c(.24,.01,.04,.01,.25,-.02,.04,-.02,.24),3,3,
                dimnames=list(names(target),names(target)))
Z <- as.matrix(dat[,cols])
wf <- get.weight(Z,target)
args <- list(data=dat,w=wf$weight,model_covariates="z2",effect_modifiers="z1",
             return_full=TRUE)
fixed <- suppressWarnings(do.call(fit_robust_lb,args))
lb <- suppressWarnings(fit_logbinomial(dat,wf$weight,model_covariates="z2",
                                      effect_modifiers="z1"))
stopifnot(identical(fixed$point.est,coef(lb)),
          identical(names(fixed$point.est),c("(Intercept)","x","z2","x:z1")),
          identical(dim(fixed$vcov),c(4L,4L)))

## 固定权重版直接使用 vcovHC；核对 HC0、HC1 和 HC3 的完整协方差。
for (hc in c("HC0","HC1","HC3")) {
  result <- suppressWarnings(do.call(fit_robust_lb,c(args,list(hc_type=hc))))
  expected <- sandwich::vcovHC(lb,type=hc)
  stopifnot(identical(result$point.est,coef(lb)),
            max(abs(result$vcov-expected))<1e-12,
            abs(result$se^2-expected["x","x"])<1e-12)
}

## 联合版的独立对照：逐人 log likelihood 的数值 score 和观测 Hessian。
X <- model.matrix(lb)
par <- coef(lb)
loglik.rows <- function(b) {
  p <- exp(drop(X%*%b))
  dat$y*log(p)+(1-dat$y)*log1p(-p)
}
numerical.score <- function(b) sapply(seq_along(b),function(j) {
  step <- rep(0,length(b)); step[j] <- 1e-6
  (loglik.rows(b+step)-loglik.rows(b-step))/(2e-6)
})
score <- wf$weight*numerical.score(par)
H <- optimHess(par,function(b) -sum(wf$weight*loglik.rows(b)),
               control=list(ndeps=rep(1e-4,length(par))))
V.observed <- solve(H)%*%crossprod(score)%*%solve(H)
scaled <- suppressWarnings(do.call(fit_robust_lb,modifyList(args,list(w=7*wf$weight))))
stopifnot(max(abs(scaled$point.est-fixed$point.est))<1e-7,
          max(abs(scaled$vcov-fixed$vcov))<1e-7)
cat("PASS: fixed-weight RLB matches vcovHC HC0/HC1/HC3 and is invariant to weight rescaling\n")

jargs <- c(args,list(gamma=wf$ga,target_mean=target,weight_covariates=cols))
joint <- suppressWarnings(do.call(fit_robust_lb_joint,jargs))
full <- suppressWarnings(do.call(fit_robust_lb_joint,c(jargs,
  list(target_cov=Sigma,m_AD=300,account_for_target_mean=TRUE))))
off <- suppressWarnings(do.call(fit_robust_lb_joint,c(jargs,list(account_for_gamma=FALSE))))
stopifnot(identical(joint$point.est,fixed$point.est),
          identical(full$point.est,fixed$point.est),
          max(abs(off$vcov-V.observed))<1e-7)

## 独立对照：数值求整个 (gamma,beta) 估计方程的 Jacobian。
q <- ncol(Z); n <- nrow(dat); p <- ncol(X)
stacked <- function(a) {
  D <- sweep(Z,2,target[cols],"-")
  w <- exp(drop(D%*%a[seq_len(q)])); w <- w/mean(w)
  cbind(w*D,w*numerical.score(a[q+seq_len(p)]))
}
a <- c(wf$ga,par)
U <- stacked(a)
A <- sapply(seq_along(a),function(j) {
  step <- rep(0,length(a)); step[j] <- 1e-4
  (colSums(stacked(a+step))-colSums(stacked(a-step)))/(2e-4)
})
Ai <- solve(A)
V.joint <- Ai%*%crossprod(U)%*%t(Ai)
stopifnot(max(abs(joint$vcov-V.joint[q+seq_len(p),q+seq_len(p)]))<1e-7)
cat("PASS: joint RLB matches independently differentiated stacked estimating equations\n")

## AD 协方差转换、按名称重排，以及样本量缩放。
twice <- suppressWarnings(do.call(fit_robust_lb_joint,c(jargs,
  list(target_cov=Sigma,m_AD=600,account_for_target_mean=TRUE))))
direct <- suppressWarnings(do.call(fit_robust_lb_joint,c(jargs,
  list(target_mean_vcov=Sigma[c(3,2,1),c(3,2,1)]/300,account_for_target_mean=TRUE))))
delta <- full$vcov-joint$vcov
stopifnot(max(abs(twice$vcov-joint$vcov-delta/2))<1e-12,
          max(abs(direct$vcov-full$vcov))<1e-12,
          max(abs(delta))>1e-6,
          min(eigen(delta,symmetric=TRUE)$values)>-1e-10,
          max(abs(full$vcov-full$vcov.ipd-full$vcov.target))<1e-12)

## 扰动 AD 均值并重拟合，用独立有限差分核对 AD 的方差贡献。
refit <- function(mu) {
  w <- apply(sweep(Z,2,mu,function(z,u) ifelse(z==1,2*u,2*(1-u))),1,prod)
  suppressWarnings(coef(glm(y~x+z2+x:z1,data=dat,family=binomial(link="log"),
    weights=w,start=par,control=glm.control(epsilon=1e-12,maxit=100))))
}
h <- 5e-4
J <- sapply(seq_len(q),function(j) {
  step <- rep(0,q); step[j] <- h
  (refit(target[cols]+step)-refit(target[cols]-step))/(2*h)
})
V.target <- J%*%(Sigma[cols,cols]/300)%*%t(J)
stopifnot(max(abs(V.target-full$vcov.target))/max(abs(V.target))<.001)
cat("PASS: AD covariance propagation matches target-mean perturbation/refits and 1/m scaling\n")

## 无平衡列、无协变量、分数列、多个 modifier 及单效应返回接口。
empty <- suppressWarnings(fit_robust_lb_joint(dat,weight_covariates=character(0),
  model_covariates=character(0),account_for_target_mean=TRUE,return_full=TRUE))
simple <- suppressWarnings(fit_robust_lb(dat,model_covariates=character(0),return_full=TRUE))
stopifnot(identical(empty$point.est,simple$point.est),
          max(abs(empty$vcov-simple$vcov))<1e-7)
arm.var <- sapply(0:1,function(x) {
  y <- dat$y[dat$x==x]; mu <- mean(y)
  sum((y-mu)^2)/(length(y)^2*mu^2)
})
stopifnot(abs(simple$se^2-sum(arm.var))<1e-9)
dat$h_E <- (dat$z1-.4)/sqrt(.24)
score.fit <- suppressWarnings(fit_robust_lb_joint(dat,
  weight_covariates=character(0),model_covariates=c("h_E","z2"),
  effect_modifiers=c("h_E","z2"),return_full=TRUE))
stopifnot(length(score.fit$point.est)==6L,all(is.finite(score.fit$vcov)))
cat("PASS: zero weight columns, intercept-only, score columns and multiple modifiers\n")

## 比较入口的矩阵/完整列表都包含新增方法；失败占位不会丢掉其他方法。
ipd <- list(data=dat)
ad <- get.aggre(ipd); ad$anchored.mean.x <- target; ad$anchored.cov.x <- Sigma
comparison <- suppressWarnings(get.compare.anchored(ipd,ad,run_exact=FALSE,
  weight_covariates=cols,model_covariates="z2",effect_modifiers="z1"))
stopifnot("RLB" %in% names(comparison$fits),
          "brm-FC" %in% names(comparison$fits),
          sum(comparison$estimates$method=="RLB")==4L)
scalar <- suppressWarnings(get.compare.anchored(ipd,ad,run_exact=FALSE,
  weight_covariates=cols,model_covariates="z2"))
stopifnot(identical(dim(scalar),c(5L,length(anchored_method_names()))),
          all(c("brm-FC","RLB",anchored_glm_fc_method_names()) %in% colnames(scalar)))
failed <- comparison_parameter_rows("RLB",NULL,"robust_glm")
stopifnot(failed$status=="fit_failed",is.na(failed$estimate))
cat("PASS: comparison retains joint RLB, replaces fixed RLB with FC Wald\n")
