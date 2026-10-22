MLEst = function(group,param, y, x, va, vb, weight,gamma.weight, max.step, thres, alpha.start,
    beta.start, pa, pb,target_mean, method = "brm", CI = "wald",
    target_var = NULL, m_AD = NULL, target_cov = NULL,
    target_mean_vcov = NULL, vg = NULL, optim.maxit = 200L, optim.reltol = 1e-6) {

    ## starting values for parameter optimization
    if (is.null(alpha.start))
        alpha.start = rep(0, pa)
    if (is.null(beta.start))
        beta.start = rep(0, pb)

    if (param == "OR") {
        fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                         weight = weight, start = c(beta.start, alpha.start))

        point.temp = summary(fit)$coefficients[, 1]
        index = c((pb + 1):(pa + pb), 1:pb)
        point.est = point.temp[index]

        cov = stats::vcov(fit)[index, index]

        converged = fit$converged

    } else {

      ### point estimate
      if(method == "brm"){
          mle = max.likelihood(param, y, x, va, vb, alpha.start, beta.start,
              weight, max.step, thres, pa, pb, optim.maxit = optim.maxit,
              optim.reltol = optim.reltol)
      } else if(method == "firth"){
          if(param == "RR"){
              mle = max.likelihood.firth.rr(param, y, x, va, vb, alpha.start, beta.start,
                  weight, max.step, thres, pa, pb)
          } else {
              mle = max.likelihood.firth.rd(param, y, x, va, vb, alpha.start, beta.start,
                  weight, max.step, thres, pa, pb)
          }
      } else if(method == "jeffrey-p"){
            mle = max.likelihood.jeffrey.direct(param, y, x, va, vb, alpha.start, beta.start,
                  weight, max.step, thres, pa, pb)
      } else if(method == "jeffrey-est"){
            mle = max.likelihood.jeffrey(param, y, x, va, vb, alpha.start, beta.start,
                  weight, max.step, thres, pa, pb)
      } else {
          stop(paste0(method, " is not a recognized method!"))
      }

      point.est = mle$par
      converged = mle$convergence
      if (!isTRUE(converged)) {
        name <- paste(c(rep("alpha",pa),rep("beta",pb)),c(seq_len(pa),seq_len(pb)))
        candidate <- WrapResults(point.est,matrix(NA_real_,pa+pb,pa+pb),param,name,va,vb,FALSE)
        candidate$last_optimizer <- mle
        candidate$optimization <- mle$diagnostics
        candidate$method <- method
        candidate$optimization_input <- list(param=param,y=y,x=x,va=va,vb=vb,weight=weight,
          alpha.start=point.est[seq_len(pa)],beta.start=point.est[pa+seq_len(pb)],
          max.step=max.step,thres=thres,pa=pa,pb=pb)
        if (method == "brm") {
          candidate$optimization_input$optim.maxit <- optim.maxit
          candidate$optimization_input$optim.reltol <- optim.reltol
        }
        if (method == "firth") candidate$fc.diagnostics <- mle$diagnostics
        warning(method," did not converge; last iterate retained separately for inspection/retry.",call.=FALSE)
        return(candidate)
      }
      alpha.ml = point.est[1:pa]
      beta.ml = point.est[(pa + 1):(pa + pb)]
      if(group=="AD"){
        if (method == "firth") {
          derivatives = if (param == "RR") firth_rr_derivatives else firth_rd_derivatives
          terms = derivatives(alpha.ml, beta.ml, x, va, vb)
          ## Fixed-weight, model-based FIRST-ORDER sandwich; not inverse Fisher.
          cov = firth_score_adjustment(y, va, vb, weight, terms)$covariance
        } else if (param == "RR") {
          cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb,weight)
        } else if (param == "RD") {
          cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb,weight)
        }
        sd.est = sqrt(diag(cov))
      }else if(group=="IPD"){
        gamma.ml = gamma.weight
        if (is.null(vg)) vg = vb[, -1, drop = FALSE]
        if (method == "firth") {
          ## The joint sandwich reconstructs weights from gamma. It must use
          ## the same relative weights as the FC point estimate.
          joint.weight = normalize_exp_weights(sweep(vg, 2, target_mean, "-") %*% gamma.ml)
          fit.weight = weight / max(weight)
          fit.weight = fit.weight / mean(fit.weight)
          if (length(joint.weight) != length(weight) || any(!is.finite(joint.weight)) ||
              max(abs(joint.weight - fit.weight)) > 1e-6 * max(1, fit.weight)) {
            stop("FC joint variance requires gamma.weight/gamma.ipd matching the supplied weights; use the fitted MAIC gamma, or group='AD' for fixed-weight covariance.")
          }
        }
        if (param == "RR") {
          cov.joint = tryCatch(
            var.mle.rr.joint(
              y,x, alpha.ml, beta.ml,gamma.ml, va, vb,vg, target_mean,
              target_var = target_var,
              m_AD = m_AD,
              target_cov = target_cov,
              target_mean_vcov = target_mean_vcov
            ),
            error = function(e) {
              warning("Joint sandwich with AD mean uncertainty failed: ",
                      conditionMessage(e),
                      "; retrying fixed-target joint sandwich.",
                      call. = FALSE)
              var.mle.rr.joint(y,x, alpha.ml, beta.ml,gamma.ml,
                                va, vb,vg, target_mean)
            }
          )
        }
        if (param == "RD") {
          cov.joint = tryCatch(
            var.mle.rd.joint(
              y,x, alpha.ml, beta.ml,gamma.ml, va, vb,vg, target_mean,
              target_var = target_var,
              m_AD = m_AD,
              target_cov = target_cov,
              target_mean_vcov = target_mean_vcov
            ),
            error = function(e) {
              warning("Joint sandwich with AD mean uncertainty failed: ",
                      conditionMessage(e),
                      "; retrying fixed-target joint sandwich.",
                      call. = FALSE)
              var.mle.rd.joint(y,x, alpha.ml, beta.ml,gamma.ml,
                                va, vb,vg, target_mean)
            }
          )
        }
        q = length(gamma.ml)
        idx = (q + 1):(q + pa + pb)
        cov = cov.joint[idx, idx, drop = FALSE]
        if (any(!is.finite(diag(cov))) || any(diag(cov) < 0)) {
          warning("Joint sandwich produced non-finite or negative variances; retrying fixed-target joint sandwich.",
                  call. = FALSE)
          if (param == "RR") {
            cov.joint = var.mle.rr.joint(y,x, alpha.ml, beta.ml,gamma.ml,
                                         va, vb,vg, target_mean)
          }
          if (param == "RD") {
            cov.joint = var.mle.rd.joint(y,x, alpha.ml, beta.ml,gamma.ml,
                                         va, vb,vg, target_mean)
          }
          cov = cov.joint[idx, idx, drop = FALSE]
        }
        sd.est = sqrt(diag(cov))
      }
      

      ### Computing Fisher Information:
      
    }

    conf.lower = point.est + stats::qnorm(0.025) * sd.est
    conf.upper = point.est + stats::qnorm(0.975) * sd.est
    p.temp = stats::pnorm(point.est/sd.est, 0, 1)
    p.value = 2 * pmin(p.temp, 1 - p.temp)

    if(CI == "wald"){
      ci.est = list(low = conf.lower,
                    up = conf.upper,
                    p = p.value)
    }else if(CI == "exact"){
      ci.est.alpha = exact(param,y, x, va, vb, weight, max.step, thres, thres.dicho=1e-3, point.est, sd.est, pa, pb)
      ci.est = list(low = c(ci.est.alpha$low,conf.lower[(pa+1):(pa+pb)]),
                    up = c(ci.est.alpha$up,conf.upper[(pa+1):(pa+pb)]),
                    p = c(ci.est.alpha$p,p.value[(pa+1):(pa+pb)]))
    }else if(CI == "LRT"){
      ci.est.alpha = profile(param,y, x, va, vb, weight, max.step, thres, point.est, sd.est, pa, pb)
      ci.est = list(low = c(ci.est.alpha$low,conf.lower[(pa+1):(pa+pb)]),
                    up = c(ci.est.alpha$up,conf.upper[(pa+1):(pa+pb)]),
                    p = c(ci.est.alpha$p,p.value[(pa+1):(pa+pb)]))
    }else {
      stop(paste0(CI, " is not a recognized method!"))
    }

    name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
    sol = WrapResults(point.est, cov, param, name, va, vb, converged, ci.est)
    if (method == "brm" && param %in% c("RR", "RD")) sol$optimization = mle$diagnostics
    if (method == "firth" && param %in% c("RR", "RD")) {
      sol$fc.diagnostics = mle$diagnostics
      sol$fc.weight.interpretation = "fixed_analytic"
      sol$variance.method = if (group == "AD") "fixed_weight_model_sandwich_first_order" else "joint_sandwich_first_order"
    }
    return(sol)
}

