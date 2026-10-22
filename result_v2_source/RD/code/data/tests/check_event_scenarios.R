## 在 anchored_v1 运行：Rscript data/tests/check_event_scenarios.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored_data()

common <- joint_rr_gop_design(event = "Common", ess_ratio = .50)
rare <- joint_rr_gop_design(event = "rare", ess_ratio = .50)

## 核对用户约定的总体风险；容差只允许表格的四位小数舍入误差。
targets <- list(common = c(.4010, .5990, .4186, .5022, 1.1997),
                rare = c(.0501, .0763, .0488, .0589, 1.2070))
for (d in list(common, rare)) {
  t <- d$marginal_truth
  actual <- c(t["AC", "p_C"], t["AC", "p_A"],
              t["BC", "p_C"], t["BC", "p_B"], exp(t["BC", "log_RR_BC"]))
  stopifnot(max(abs(actual - targets[[d$event]])) < 5e-5)
  set.seed(23)
  sim <- generate_anchored_data(100, 150, d)
  stopifnot(identical(sim$design$event, d$event),
            all(sim$data.IPD$count[c("Na0", "Na1")] == 50),
            all(sim$data.AD.or$count[c("Na0", "Na1")] == 75))
}
stopifnot(identical(common$risk_table$log_RR_AC, rare$risk_table$log_RR_AC),
          identical(common$risk_table$log_RR_BC, rare$risk_table$log_RR_BC),
          all(rare$risk_table$p_C < common$risk_table$p_C),
          all(rare$risk_table$p_A < common$risk_table$p_A),
          all(rare$risk_table$p_B < common$risk_table$p_B))
cat("PASS: common/rare marginal targets, unchanged conditional RRs and sampling\n")

## 原位置调用给出同一风险表；场景选择忽略大小写，显式参数优先。
manual <- joint_rr_gop_design(.20, .19, .30, .50, -.26, .20, .40, ess_ratio=.50)
stopifnot(identical(common$risk_table, manual$risk_table),
          identical(common$marginal_truth, manual$marginal_truth),
          identical(common, joint_rr_gop_design(event="COMMON", ess_ratio=.50)),
          identical(rare, joint_rr_gop_design(event="Rare", ess_ratio=.50)))
override <- joint_rr_gop_design(event="rare", alpha_AC=.40, rho=0, beta0=-3,
                                ess_ratio=.50)
reference <- joint_rr_gop_design(.40, .19, .30, 0, -3, .20, .40, ess_ratio=.50)
stopifnot(identical(override$risk_table, reference$risk_table),
          max(abs(override$marginal_truth[, "log_RR_BC"] - .19)) < 1e-10)
q <- c(.2, .4, .6, .8, .5, .3)
direct <- joint_rr_gop_design(event="Common", q=q)
stopifnot(identical(unname(direct$q), q),
          identical(direct$risk_table, common$risk_table),
          joint_rr_gop_design(event="Common")$ess_ratio == 1)
cat("PASS: manual coefficients, case handling, explicit overrides and population controls\n")

## 拼写错误及不完整的手动参数不得静默采用其他场景。
for (bad in list("unknown", "com", "", NA_character_, character(0),
                 c("common", "rare"), 1)) {
  stopifnot(inherits(try(joint_rr_gop_design(event=bad), silent=TRUE), "try-error"))
}
stopifnot(inherits(try(joint_rr_gop_design(alpha_AC=.20), silent=TRUE), "try-error"))
cat("PASS: invalid scenarios and missing coefficients are rejected\n")
