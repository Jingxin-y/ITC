## 联合 RR-GOP 设计及共享抽样器：先固定设计和公共风险表，再抽样两项试验。
## 共享抽样器也支持 webster_rr_design() 和独立的 webster_rd_design()；
## 三种机制使用同样的分数和 16 行风险表接口。估计流程仍按 RR 运行。
## 在 anchored_v1 目录运行：
## source("load_anchored.R", encoding="UTF-8"); source_anchored_data()
## 依赖 data/helpers/joint_rr_gop.R（ESS、分数、风险表、真值积分）、
## data/helpers/getProbRR_GOP.R（三臂概率）、data/ad_summary.R（AD 汇总）及 data/helpers/covariates.R。
## 这些函数按模块分开定义，不是仅 source 本文件就会加载。
## event 选择两位小数的 common/rare 参数；不区分大小写。
joint_rr_gop_parameters <- function(event) {
  if (!is.character(event) || length(event) != 1L || is.na(event) ||
      !tolower(event) %in% c("common", "rare")) {
    stop('event must be "common" or "rare" (case-insensitive).')
  }
  list(alpha_AC = 0.20, alpha_BC = 0.19, lambda = 0.30, rho = 0.50,
       beta0 = if (tolower(event) == "common") -0.26 else -8.53,
       beta_E = 0.20, beta_P = 0.40)
}

## 指定 event 时，未提供的系数取场景值；显式系数优先，可用于敏感性分析。
## 不指定 event 时，仍须给出全部七个系数，兼容原来的具名和位置调用。
## 约 0.40/0.60 与 0.05/0.076 的 AC 风险对应 ess_ratio=0.50 和默认 direction。
joint_rr_gop_design <- function(alpha_AC = NULL, alpha_BC = NULL,
                                lambda = NULL, rho = NULL,
                                beta0 = NULL, beta_E = NULL, beta_P = NULL,
                                q = NULL, ess_ratio = 1,
                                direction = c(-1, 1, -1, 1, -1, 1),
                                event = NULL) {
  parameters <- list(alpha_AC = alpha_AC, alpha_BC = alpha_BC,
                     lambda = lambda, rho = rho,
                     beta0 = beta0, beta_E = beta_E, beta_P = beta_P)
  if (!is.null(event)) {
    preset <- joint_rr_gop_parameters(event)
    supplied <- !vapply(parameters, is.null, logical(1))
    preset[supplied] <- parameters[supplied]
    parameters <- preset
    event <- tolower(event)
  }
  missing_parameters <- names(parameters)[vapply(parameters, is.null, logical(1))]
  if (length(missing_parameters)) {
    stop('Specify event="common"/"rare" or provide all coefficients. Missing: ',
         paste(missing_parameters, collapse = ", "))
  }
  e <- c(.7, .3, .7, .3, .7, .3)  # 固定 AD / BC 概率
  ## 由 ESS 及给定方向求 IPD / AC 概率；直接指定 q 的入口也保留。
  if (is.null(q)) {
    q <- ipd_prob_for_ess_theory(ess_ratio = ess_ratio, e = e, direction = direction)
  }
  stopifnot(length(q) == 6L, all(q >= 0 & q <= 1), parameters$rho %in% c(0, .5, 1))
  names(e) <- names(q) <- paste0("z", 1:6)
  risk_table <- do.call(rr_gop_risk_table, parameters)

  list(parameters = parameters, e = e, q = q,
       ess_ratio = rel_ess_binary(q, e),
       score_center = .9, score_scale = sqrt(.63),
       risk_table = risk_table,
       marginal_truth = rbind(AC = rr_gop_marginal_truth(risk_table, q),
                             BC = rr_gop_marginal_truth(risk_table, e)),
       event = event)
}

## group="IPD" 表示 A-C 试验；group="AD" 表示 B-C 试验。
## x=0 为 C，x=1 为该试验的活性治疗；保留原 data/count/design 返回结构。
## n 为正偶数；随机分配后两组各有 n/2 人。
data.generation <- function(group, n, design) {
  group <- match.arg(group, c("IPD", "AD"))
  if (n < 2 || n %% 2 != 0) stop("n must be a positive even integer for equal allocation.")
  z.prob <- if (group == "IPD") design$q else design$e  # q 变化，e 始终固定
  z <- binary_covariates(n, z.prob)
  scores <- rr_gop_scores(z)
  index <- 1L + scores$S_E + 4L * scores$S_P
  table <- design$risk_table
  p0.true <- table$p_C[index]
  p1.true <- if (group == "IPD") table$p_A[index] else table$p_B[index]

  x <- sample(rep(0:1, each = n / 2))
  pA.true <- p0.true
  pA.true[x == 1] <- p1.true[x == 1]
  y <- rbinom(n, 1, pA.true)
  count <- c(Na0 = sum(x == 0), Na1 = sum(x == 1),
             N0_1 = sum(y[x == 0]), N1_1 = sum(y[x == 1]))

  list(
    data = data.frame(y = y, x = x, z, scores, check.names = FALSE),
    count = count,
    design = list(trial = if (group == "IPD") "AC" else "BC",
                  target_mean = design$e, ess_ratio = design$ess_ratio,
                  n_cov = 6L, covariate_prob = z.prob),
    truth = list(p0 = p0.true, p1 = p1.true, observed_risk = pA.true,
                 risk_table_row = index)
  )
}

## 同一 design 传给两项试验；不同重复也可继续使用这张风险表。
generate_anchored_data <- function(n_AC, n_BC = n_AC, design) {
  data.IPD <- data.generation("IPD", n_AC, design)
  data.AD.or <- data.generation("AD", n_BC, design)
  list(data.IPD = data.IPD, data.AD.or = data.AD.or,
       data.AD = get.aggre(data.AD.or), design = design)
}
