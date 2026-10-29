## 三臂联合风险转换；输入均为对数尺度，GOP 是 A/B/C 三臂 odds 的乘积。
## 给定 log(RR_AC)、log(RR_BC)、log(GOP_ABC)，返回 p_C、p_A、p_B。
## 此处保留原 rr_gop_risk_table() 的求根代码，独立供数据生成调用。
getProbScalarRR_GOP <- function(logrr_AC, logrr_BC, loggop_ABC) {
  log_R <- c(0, logrr_AC, logrr_BC)
  offset <- log_R - max(log_R)

  ## c = exp(-max(log_R))*expit(u)，保证 0<c<1/max(1,RA,RB)。
  ## 方程为 logit(p_C)+logit(p_A)+logit(p_B)=log(GOP_ABC)。
  equation <- function(u) {
    log_p <- offset + plogis(u, log.p = TRUE)
    sum(log_p - log(-expm1(log_p))) - loggop_ABC
  }
  root <- uniroot(equation, c(-30, 30), extendInt = "upX", tol = 1e-12)$root
  setNames(exp(offset + plogis(root, log.p = TRUE)), c("p_C", "p_A", "p_B"))
}

## 向量版本：三个输入长度相同，或其中某个为固定标量。
## 即使只输入一个组合，也返回具有三列的矩阵。
getProbRR_GOP <- function(logrr_AC, logrr_BC, loggop_ABC) {
  predictors <- list(logrr_AC, logrr_BC, loggop_ABC)
  n <- max(lengths(predictors))
  stopifnot(all(lengths(predictors) %in% c(1L, n)))
  predictors <- lapply(predictors, function(x) rep_len(as.numeric(x), n))
  probs <- t(vapply(seq_len(n), function(i) {
    getProbScalarRR_GOP(predictors[[1]][i], predictors[[2]][i], predictors[[3]][i])
  }, numeric(3)))
  colnames(probs) <- c("p_C", "p_A", "p_B")
  probs
}
