## 在 anchored_v1 目录运行。RD 的独立例子；原来的批量 RR 场景不变。
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
set.seed(12)

dat <- expand.grid(id = 1:100, x = 0:1, z1 = 0:1, z2 = 0:1)
v <- cbind(1, dat$z1, dat$z2)
## atanh RD(z)=0.1+0.15*z1-0.1*z2；log OP(z)=-1+0.4*z1-0.2*z2。
alpha.true <- c(.1, .15, -.1)
probs <- getProbRD(v %*% alpha.true, v %*% c(-1, .4, -.2))
dat$y <- rbinom(nrow(dat), 1, probs[cbind(seq_len(nrow(dat)),dat$x+1)])
data.IPD <- list(data=dat)
target.mean <- c(z1=.4,z2=.6)
weight.fit <- get.weight(covariate_matrix(data.IPD$data),target.mean)

fit.rd <- get.estimate(
  "IPD",data.IPD,weight.fit$weight,weight.fit$ga,target_mean=target.mean,
  method="firth",effect_modifiers=c("z1","z2"),param="RD"
)
pa <- ncol(fit.rd$va)
print(fit.rd$coefficients[seq_len(pa),,drop=FALSE]) # alpha 位于 atanh(RD) 尺度
print(fit.rd$converged)
print(head(fit.rd$param.est)) # 每个协变量取值处的条件 RD，已经取 tanh

## 可选 BC：逐个 alpha 的区间，并非整个目标人群平均 RD 的区间。
run_bc <- FALSE
if (run_bc) {
  set.seed(1000012)
  ci.rd <- exact_safe(
    "RD",dat$y,dat$x,fit.rd$va,fit.rd$vb,weight.fit$weight,
    fit.rd$point.est,fit.rd$se.est,pa,ncol(fit.rd$vb),label="RD Firth BC",
    max.step=min(pa*20,1000)
  )
  print(cbind(low=ci.rd$low,up=ci.rd$up,p=ci.rd$p))
  ## 只有 alpha 截距可这样转换为 z1=z2=0 时的 RD 区间。
  print(tanh(c(ci.rd$low[1],ci.rd$up[1])))
}
