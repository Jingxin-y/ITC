## 联合 RR-GOP 的辅助计算；z1,...,z6 对应文中的 X1,...,X6。
## 二元协变量与相对 ESS 公式沿用原 generate_RR.R。
binary_covariates <- function(n, prob) {
  z <- matrix(sapply(prob, function(p) rbinom(n, 1, p)),
              nrow = n, ncol = length(prob))
  colnames(z) <- paste0("z", seq_len(ncol(z)))
  z
}

## p 为来源 AC/IPD 分布，q 为目标 BC 分布；两者顺序不可交换。
rel_ess_binary <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  1 / prod(q^2 / p + (1 - q)^2 / (1 - p))
}

## e 是固定的 AD/BC 概率；沿 direction 移动 IPD/AC 的 q 来满足指定 ESS。
ipd_prob_for_ess_theory <- function(ess_ratio,
                                  e = c(.7, .3, .7, .3, .7, .3),
                                  direction = c(-1, 1, -1, 1, -1, 1)) {
  stopifnot(ess_ratio > 0, ess_ratio <= 1, length(direction) == 6L,
            any(direction != 0))
  q_of_s <- function(s) plogis(qlogis(e) + s * direction)
  if (ess_ratio == 1) return(e)

  f <- function(s) rel_ess_binary(q_of_s(s), e) - ess_ratio
  upper <- 1
  while (f(upper) > 0) upper <- upper * 2
  q_of_s(uniroot(f, c(0, upper), tol = 1e-12)$root)
}

rr_gop_scores <- function(z) {
  S_E <- (1 - z[, 1]) + z[, 2] + (1 - z[, 3])
  S_P <- z[, 4] + (1 - z[, 5]) + z[, 6]
  data.frame(S_E = S_E, S_P = S_P,
             h_E = (S_E - .9) / sqrt(.63),
             h_P = (S_P - .9) / sqrt(.63))
}

rr_gop_risk_table <- function(alpha_AC, alpha_BC, lambda, rho,
                              beta0, beta_E, beta_P) {
  tab <- expand.grid(S_E = 0:3, S_P = 0:3)
  tab$h_E <- (tab$S_E - .9) / sqrt(.63)
  tab$h_P <- (tab$S_P - .9) / sqrt(.63)
  tab$log_RR_AC <- alpha_AC + lambda * tab$h_E
  tab$log_RR_BC <- alpha_BC + rho * lambda * tab$h_E
  tab$log_GOP <- beta0 + beta_E * tab$h_E + beta_P * tab$h_P

  ## 三臂联合转换；每个分数组合只求一次根，两项试验随后查同一张表。
  probs <- getProbRR_GOP(tab$log_RR_AC, tab$log_RR_BC, tab$log_GOP)
  tab$p_C <- probs[, "p_C"]
  tab$p_A <- probs[, "p_A"]
  tab$p_B <- probs[, "p_B"]
  tab
}

## 独立 Bernoulli 之和的概率；AC 的三项成功概率未必相等。
binary_score_prob <- function(prob) {
  mass <- 1
  for (p in prob) mass <- c(mass * (1 - p), 0) + c(0, mass * p)
  mass
}

rr_gop_marginal_truth <- function(risk_table, covariate_prob) {
  e_prob <- binary_score_prob(c(1 - covariate_prob[1], covariate_prob[2],
                                1 - covariate_prob[3]))
  p_prob <- binary_score_prob(c(covariate_prob[4], 1 - covariate_prob[5],
                                covariate_prob[6]))
  mass <- e_prob[risk_table$S_E + 1] * p_prob[risk_table$S_P + 1]
  risk <- colSums(risk_table[c("p_A", "p_B", "p_C")] * mass)
  c(risk, log_RR_AC = unname(log(risk["p_A"] / risk["p_C"])),
    log_RR_BC = unname(log(risk["p_B"] / risk["p_C"])),
    log_RR_AB = unname(log(risk["p_A"] / risk["p_B"])))
}
