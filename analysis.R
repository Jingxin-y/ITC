rm(list=ls())
gc()

setwd("D:/UofT/code/RRRDOR/ITC/result/anchored")

library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(matrixStats)
library(stringr)

est.result <- function(df, para.true) {
  est <- df$estimate
  bias <- mean(est) - para.true
  
  se <- df$se
  se.est <- mean(se) / sqrt(length(se))
  
  sd.est <- mean(se)
  sd.mc <- sd(est)
  acc <- sd.est / sd.mc
  # 
  low <- df$low
  up <- df$up
  # low <- est - 8*se
  # up <- est + 8*se
  cov <- mean((low < para.true) * (up > para.true))
  
  p <- mean(df$p <= 0.05)
  
  estm.ml <- c(bias, se.est, acc, cov, p)
  return(estm.ml)
}

param = 'RR'
event = 'common'
hypothesis = 'alternative'
n=50
ess = 0.35


data1 <- read.csv(paste0("ITC_simulation_results_anchored_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 200,"_ess_",ess,".csv"))
data2 <- read.csv(paste0("ITC_simulation_results_anchored_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 400,"_ess_",ess,".csv"))
data3 <- read.csv(paste0("ITC_simulation_results_anchored_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 600,"_ess_",ess,".csv"))
data4 <- read.csv(paste0("ITC_simulation_results_anchored_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 800,"_ess_",ess,".csv"))
data5 <- read.csv(paste0("ITC_simulation_results_anchored_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 1000,"_ess_",ess,".csv"))
data = rbind(data1,data2,data3,data4,data5)

{
  data.brm <- data[(seq_len(nrow(data)) %% 11) == 1,]
  data.brm.an <- data[(seq_len(nrow(data)) %% 11) == 2,]
  data.lb.an <- data[(seq_len(nrow(data)) %% 11) == 3,]
  data.lp.an <- data[(seq_len(nrow(data)) %% 11) == 4,]
  data.rlp.an <- data[(seq_len(nrow(data)) %% 11) == 5,]
  data.CMH.an <- data[(seq_len(nrow(data)) %% 11) == 6,]
  data.RCMH <- data[(seq_len(nrow(data)) %% 11) == 7,]
  data.brm.ad.an <- data[(seq_len(nrow(data)) %% 11) == 8,]
  data.GC <- data[(seq_len(nrow(data)) %% 11) == 9,]
  data.exact.an <- data[(seq_len(nrow(data)) %% 11) == 10,]
  data.exact.ad.an <- data[(seq_len(nrow(data)) %% 11) == 0,]
  
  data.RLP <- read.csv(paste0("ITC_simulation_results_anchored_RLP_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", 1000,"_ess_",ess,".csv"))
  data.rlp.an <- data.RLP[(seq_len(nrow(data.RLP)) %% 11) == 5,]
  num_cols <- sapply(data.lb.an, is.numeric)
  df.lb.an <- data.lb.an[ apply(data.lb.an[ , num_cols], 1, function(x) all(is.finite(x))), ]
  num_cols <- sapply(data.lp.an, is.numeric)
  df.lp.an <- data.lp.an[ apply(data.lp.an[ , num_cols], 1, function(x) all(is.finite(x))), ]
  
}


if(param == 'RR'){
  if(event == 'common'){          
    if(hypothesis == 'alternative'){
      alpha.true <- 0.4         ### 0.5,0.75
      beta.true  <- c(1.3, -0.5)
      gamma.true <- c(0,0)
    }else{
      alpha.true <- 0         ### 0.7
      beta.true <- c(1.5, 0.6)
      gamma.true <- c(0, 0)
    }
  }else{                         
    if(hypothesis == 'alternative'){
      alpha.true <- 0.7     ### 0.05,0.1
      beta.true  <- c(-5, -0.5)
      gamma.true <- c(0,0)
    }else{
      alpha.true <- 0          ### 0.09
      beta.true <- c(-4.7, 0.5)
      gamma.true <- c(0, 0)
    }
  }
}


{
  result.data.brm <- est.result(data.brm,alpha.true)
  result.data.brm.an <- est.result(data.brm.an,alpha.true)
  result.data.lb.an <- est.result(df.lb.an,alpha.true)
  result.data.lp.an <- est.result(df.lp.an,alpha.true)
  result.data.rlp.an <- est.result(data.rlp.an,alpha.true)
  result.data.CMH.an <- est.result(data.CMH.an,alpha.true)
  result.data.RCMH <- est.result(data.RCMH,alpha.true)
  result.data.brm.ad.an <- est.result(data.brm.ad.an,alpha.true)
  result.data.GC <- est.result(data.GC,alpha.true)
  result.data.exact.an <- est.result(data.exact.an,alpha.true)
  result.data.exact.ad.an <- est.result(data.exact.ad.an,alpha.true)
  
  result <- cbind(result.data.brm,result.data.CMH.an,
                  result.data.RCMH,result.data.lb.an,result.data.lp.an,
                  result.data.rlp.an,result.data.GC,result.data.brm.an,
                  result.data.brm.ad.an,result.data.exact.an,
                  result.data.exact.ad.an)
  
  rownames(result) <- c("bias", "se", "acc", "coverage", "p")
  colnames(result) <- c("brm","CMH.an","RCMH","lb.an",
                        "lp.an","rlp.an","GC","brm.an","brm.ad.an",
                        "brm.bc.an","brm.adbc.an")
  
}

result

data_names <- c("data.brm","data.CMH.an","data.RCMH","data.lb.an",
                "data.lp.an","data.rlp.an","data.GC","data.brm.an",
                "data.brm.ad.an","data.exact.an","data.exact.ad.an")

data_list <- mget(data_names)

est <- do.call(cbind, lapply(data_list, `[[`, "estimate"))
colnames(est) <- c("brm","CMH.an","RCMH","lb.an","lp.an","rlp.an",
                   "GC","brm.an","brm.ad.an","brm.bc.an","brm.adbc.an")

se <- do.call(cbind, lapply(data_list, `[[`, "se"))
colnames(se) <-c("brm","CMH.an","RCMH","lb.an","lp.an","rlp.an",
                 "GC","brm.an","brm.ad.an","brm.bc.an","brm.adbc.an")


est <- est[, c("brm","CMH.an","RCMH","lb.an","lp.an","rlp.an",
               "GC","brm.an","brm.ad.an","brm.bc.an","brm.adbc.an"), drop = FALSE]
est_df <- as.data.frame(est) 
est_long <- melt(est, variable.name = "method", value.name = "estimate")
colnames(est_long) <- c("number","method","estimate")

dynamic_ylim <- function(x, lower = -Inf, upper = Inf, include = NULL) {
  finite <- x[is.finite(x)]
  if (length(finite) == 0) return(NULL)
  
  clip_low <- is.finite(lower) && any(finite < lower)
  clip_high <- is.finite(upper) && any(finite > upper)
  if (!clip_low && !clip_high) return(NULL)
  
  visible <- finite
  if (clip_low) visible <- visible[visible >= lower]
  if (clip_high) visible <- visible[visible <= upper]
  visible <- c(visible, include)
  visible <- visible[is.finite(visible)]
  
  low <- if (clip_low) lower else min(visible)
  high <- if (clip_high) upper else max(visible)
  if (!is.finite(low) || !is.finite(high) || low >= high) {
    if (is.finite(lower) && is.finite(upper) && lower < upper) {
      return(c(lower, upper))
    }
    return(range(finite))
  }
  
  pad <- diff(c(low, high)) * 0.06
  if (!clip_low) low <- low - pad
  if (!clip_high) high <- high + pad
  c(low, high)
}


p1 = ggplot(est_long, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(color = "grey20",width = 0.8, outlier.size = 0.5, alpha = 0.9) +
  ggsci::scale_fill_d3(palette = "category20") +
  theme_minimal(base_size = 18) +
  geom_hline(yintercept = alpha.true, linetype = "dashed", color = "steelblue") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "white", color = "black") +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +
  labs(title = paste0(
    "Monte Carlo distributions of estimated log(RR)s \nacross methods (n = ",n, ")"),
    x = "Method", y = "Estimated log(RR)")

est_ylim <- dynamic_ylim(est_long$estimate, lower = -4, upper = 4, include = alpha.true)
if (!is.null(est_ylim)) {
  p1 <- p1 + coord_cartesian(ylim = est_ylim)
}

removed <- est_long %>%
  group_by(method) %>%
  summarise(
    n_hi = sum(!is.na(estimate) & estimate >  4),   # >5
    n_lo = sum(!is.na(estimate) & estimate < -4),   # <-5
    .groups = "drop"
  )

y_hi <-  4
y_lo <- -4


p1 <- p1 +
  geom_text(
    data = filter(removed, n_hi > 0),
    aes(x = method, y = y_hi, label = n_hi),
    inherit.aes = FALSE,
    size = 4, color = "blue"
  ) +
  geom_text(
    data = filter(removed, n_lo > 0),
    aes(x = method, y = y_lo, label = n_lo),
    inherit.aes = FALSE,
    size = 4, color = "blue"
  )

p1 <- p1 +
  geom_point(
    data = filter(removed, n_hi > 0),
    aes(x = method, y = y_hi, size = n_hi),
    inherit.aes = FALSE,
    shape = 21, fill = "red", color = "red", alpha = 0.55
  ) +
  geom_point(
    data = filter(removed, n_lo > 0),
    aes(x = method, y = y_lo, size = n_lo),
    inherit.aes = FALSE,
    shape = 21, fill = "red", color = "red", alpha = 0.55
  ) +
  scale_size_area(max_size = 20)
p1
ggsave(paste0("ITC_est_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p1, width = 8, height = 6, dpi = 300)


### plot of SE
se <- se[,c("brm","CMH.an","RCMH","lb.an","lp.an","rlp.an",
            "GC","brm.an","brm.ad.an","brm.bc.an","brm.adbc.an"), drop = FALSE]
se_df <- as.data.frame(se) 
se_long <- melt(se, variable.name = "method", value.name = "SE")
colnames(se_long) <- c("number","method","SE")


se_plot <- se_long %>% filter(is.finite(SE))


removed_counts <- se_long %>%
  group_by(method) %>%
  summarise(removed = sum(!is.finite(SE) | SE > 4), .groups = "drop")

yr <- range(se_plot$SE, finite = TRUE)
y_anno <- 4
se_ylim <- dynamic_ylim(se_plot$SE, lower = 0, upper = 4)

p2 <-  ggplot(se_plot, aes(x = method, y = SE, fill = method)) +
  geom_boxplot(color = "grey20", width = 0.8, outlier.size = 0.5, alpha = 0.9) +
  ggsci::scale_fill_d3(palette = "category20") +
  theme_minimal(base_size = 18) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "white", color = "black") +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +
  labs(title = paste0(
    "Monte Carlo distributions of standard errors \nof log(RR) estimates across methods (n = ", n, ")"
  ),
  x = "Method", y = "Estimated SE")

if (!is.null(se_ylim)) {
  p2 <- p2 +
    coord_cartesian(ylim = se_ylim) +
    geom_text(
      data = dplyr::filter(removed_counts, removed > 0),
      aes(x = method, y = y_anno, label = removed),
      inherit.aes = FALSE,
      size = 4, color = "blue"
    ) +
    geom_point(
      data = dplyr::filter(removed_counts, removed > 0),
      aes(x = method, y = y_anno, size = removed),
      inherit.aes = FALSE,
      shape = 21, fill = "red", color = "red", alpha = 0.55
    ) +
    scale_size_area(max_size = 20)
}

p2
ggsave(paste0("ITC_se_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p2, width = 8, height = 6, dpi = 300)


### bar for accuracy, coverage, and p-value
method_order <- c("brm","CMH.an","RCMH","lb.an","lp.an","rlp.an",
                  "GC","brm.an","brm.ad.an","brm.bc.an","brm.adbc.an")
result <- result[, method_order, drop = FALSE]
df <- as.data.frame((result))


target <- 0.95


cov_df <- df %>%
  tibble::rownames_to_column("metric") %>%
  filter(metric == "coverage") %>%
  pivot_longer(-metric, names_to = "method", values_to = "coverage") %>%
  mutate(
    method = factor(method, levels = method_order),   
    dev = abs(coverage - target),
    label_x = ifelse(coverage >= target, coverage + 0.012, coverage - 0.012),
    hjust   = ifelse(coverage >= target, 0, 1)
  )

ylim <- range(c(cov_df$label_x, cov_df$coverage, target), finite = TRUE)
ylim <- c(ylim[1] - 0.02, ylim[2] + 0.02)

p3 <- ggplot(cov_df, aes(x = method, y = coverage)) +  
  geom_segment(
    aes(xend = method, y = target, yend = coverage, color = dev),
    linewidth = 3, lineend = "round"
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(y = label_x, label = sprintf("%.3f", coverage), hjust = hjust),
    size = 4
  ) +
  geom_hline(yintercept = target, linetype = "dashed") +
  coord_cartesian(ylim = ylim, clip = "off") +
  scale_x_discrete(limits = method_order, drop = FALSE) +  
  scale_color_gradient(low = "grey80", high = "red4") +
  guides(color = "none") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +
  labs(
    title = "Empirical coverage rates of 95% confidence intervals \nacross methods",
    x = "Method",
    y = "Coverage rate"
  )
p3
ggsave(paste0("ITC_coverage_rate_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p3, width = 8, height = 6, dpi = 300)

target <- 1

acc_df <- df %>%
  tibble::rownames_to_column("metric") %>%
  filter(metric == "acc") %>%
  pivot_longer(-metric, names_to = "method", values_to = "acc") %>%
  mutate(
    method = factor(method, levels = method_order),  # ✅固定顺序
    dev = abs(acc - target),
    label_x = ifelse(acc >= target, acc + 0.08, acc - 0.08),
    hjust   = ifelse(acc >= target, 0, 1)
  )

acc_ylim <- dynamic_ylim(acc_df$acc, lower = 0, upper = 2, include = target)
if (!is.null(acc_ylim)) {
  y_gap <- diff(acc_ylim) * 0.04
  acc_df$label_x <- pmin(pmax(acc_df$label_x, acc_ylim[1] + y_gap), acc_ylim[2] - y_gap)
  ylim <- acc_ylim
} else {
  ylim <- range(c(acc_df$label_x, acc_df$acc, target), finite = TRUE)
  ylim <- c(ylim[1] - 0.05, ylim[2] + 0.05)
}

p4 <- ggplot(acc_df, aes(x = method, y = acc)) +          # ✅不要 reorder()
  geom_segment(aes(xend = method, y = target, yend = acc, color = dev),
               linewidth = 3, lineend = "round") +
  geom_point(size = 3) +
  geom_text(aes(y = label_x, label = sprintf("%.3f", acc), hjust = hjust),
            size = 4) +
  geom_hline(yintercept = target, linetype = "dashed") +
  coord_cartesian(ylim = ylim, clip = "off") +
  scale_x_discrete(limits = method_order, drop = FALSE) + # ✅强制顺序（关键）
  scale_color_gradient(low = "grey80", high = "red4") +
  guides(color = "none") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +labs(title = "Accuracy by method", x = NULL, y = "acc")

p4
ggsave(paste0("ITC_accuracy_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p4, width = 8, height = 6, dpi = 300)


