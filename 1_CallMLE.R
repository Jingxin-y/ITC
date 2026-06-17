MLEst = function(group,param, y, x, va, vb, weight,gamma.weight, max.step, thres, alpha.start,
    beta.start, pa, pb,target_mean, method = "brm", CI = "wald",
    target_var = NULL, m_AD = NULL, target_cov = NULL,
    target_mean_vcov = NULL) {

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
              weight, max.step, thres, pa, pb)
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
      alpha.ml = point.est[1:pa]
      beta.ml = point.est[(pa + 1):(pa + pb)]
      if(group=="AD"){
        if (param == "RR") {
          cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb,weight)
        }
        if (param == "RD") {
          cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb,weight)
        }
        sd.est = sqrt(diag(cov))
      }else if(group=="IPD"){
        gamma.ml = gamma.weight
        vg = vb[, -1, drop = FALSE]
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
    return(sol)
}

