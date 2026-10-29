## Rscript estimation/tests/check_rare_boundaries.R [project]
args<-commandArgs(TRUE);project<-if(length(args))args[1] else '.'
source(file.path(project,'load_anchored.R'),encoding='UTF-8');source_anchored(project)
source(file.path(project,'estimation/cluster_all_methods.R'),encoding='UTF-8')
check<-function(ok,label){stopifnot(ok);cat('PASS:',label,'\n')}
bad<-function(expr)inherits(try(force(expr),silent=TRUE),'try-error')
v<-cbind(1,rep(c(0,1),each=20));x<-rep(0:1,20);y<-rep(c(0,0,0,1),10)
p<-c(.2,.4,-25,3)
args0<-list(y=y,x=x,va=v,vb=v,weight=rep(1,40),pars=p,va.target=v,vb.target=v,
 standardization.weight=rep(1,40),intercept=1L,alpha.bound=8,beta.bound=20,
 optim.maxit=200,optim.reltol=1e-10,score.tol=1e-4,root.tol=1e-10,
 root.maxit=200,constraint.tol=1e-7,eps=1e-12,check.deadline=function(stage)NULL)
problem<-do.call(.marginal_exact_problem,args0)
check(identical(problem$bound,c(8,20,20)) &&
 abs(problem$start[2]-(-20+1e-6))<1e-12 &&
 identical(problem$start[c(1,3)],c(.4,3)),
 'Exact repairs an out-of-bounds optimizer start while preserving other coordinates')
for(a in c(-139.3,104.1)){
 args0$pars<-c(.2,.4,a,3)
 problem<-do.call(.marginal_exact_problem,args0)
 check(all(abs(problem$start)<problem$bound) &&
  sign(problem$start[2])==sign(a) &&
  abs(problem$target(args0$pars)$value-
   marginal_logrr_target(args0$pars,v,v,args0$standardization.weight))<1e-12,
  'large signed starts are repaired without changing the original marginal target')
}
args0$pars<-c(.2,30,-5,3)
problem<-do.call(.marginal_exact_problem,args0)
check(abs(problem$start[1]-(8-1e-6))<1e-12 &&
 identical(problem$start[2:3],c(-5,3)),
 'nonintercept alpha is repaired using its explicitly supplied bound')
args0$pars<-c(.2,.4,-5,3);args0$beta.bound<-10
problem<-do.call(.marginal_exact_problem,args0)
check(identical(problem$start,c(.4,-5,3)) && identical(problem$bound,c(8,10,10)),
 'Exact retains explicitly supplied finite alpha and beta search bounds')
z<-exact_marginal_target(y=y,x=x,va=v,vb=v,pars=args0$pars,compute.interval=FALSE,
 compute.p.value=FALSE)
st<-z$diagnostics$settings
check(st$alpha.bound==15 && st$beta.bound==15 && st$optim.reltol==1e-10 &&
 st$score.tol==1e-4 && st$optim.maxit==200,
 'public Exact defaults use bounds of 15 and retain optimizer tolerances')
check(bad(validate_marginal_bc_control(list(expand.bounds=FALSE))),
 'removed expansion control is rejected by the adapter')
check(!bad(validate_marginal_bc_control(list(alpha.bound=8,beta.bound=10))),
 'pipeline adapter accepts the original bound controls')
c0<-anchored_fit_control();c1<-anchored_fit_control(list(brm_thres=1e-3,fc_thres=1e-3))
check(c0$brm_thres==1e-6 && c0$fc_thres==1e-6 && c1$brm_thres==1e-3 && c1$fc_thres==1e-3,
 'threshold sensitivity is configurable without changing defaults')
check(bad(anchored_fit_control(list(brm_thres=NA))) && bad(anchored_fit_control(list(fc_thres=0))),
 'invalid thresholds rejected')
check(c0$brm_optim_reltol==1e-6 &&
 bad(anchored_fit_control(list(brm_optim_reltol=NA))) &&
 bad(anchored_fit_control(list(brm_optim_reltol=0))),
 'inner BRM relative tolerance defaults to 1e-6 and rejects invalid values')
# The point calculation is independent of an invertible sandwich bread.
terms<-list(p0=rep(.1,40),p1=rep(.2,40),dp0=matrix(.1,40,1),dp1=matrix(.2,40,1),
 H=matrix(0,1,1),score=matrix(0,40,1),covariance=matrix(1,1,1))
point<-maic_target_point(terms,rep(1,40))
check(abs(point$estimate-log(2))<1e-12,'marginal point exists with singular bread')
zero<-maic_delta(terms,rep(1,40),matrix(rep(0:1,20),40,1),.5,matrix(.01),TRUE,rep(1,40))
check(is.finite(zero$estimate) && is.na(zero$se) && is.na(zero$low) && is.na(zero$up),
 'zero bread uses ginv but zero variance does not become a valid interval')
# The target gradient lies in the retained outcome-information subspace.
# Its influence is +/-1, orthogonal to calibration scores, so SE=1/sqrt(4).
Hnear<-matrix(c(1,1,1,1+1e-10),2,2)
signs<-c(1,-1,-1,1);zz<-c(-1,-1,1,1)
near<-list(p0=rep(.1,4),p1=rep(.2,4),dp0=matrix(0,4,2),
 dp1=matrix(.2,4,2),H=Hnear,score=signs %o% Hnear[1,],covariance=diag(2))
full<-near;full$H<-diag(c(1,2));full$score<-cbind(signs,2*signs)
out<-maic_delta(full,rep(1,4),matrix(zz,4,1),0,matrix(.01),TRUE,rep(1,4))
check(abs(out$se-1)<1e-12 && is.finite(out$low),
 'well-conditioned full-rank bread preserves the analytic SE')
out<-maic_delta(near,rep(1,4),matrix(zz,4,1),0,matrix(.01),TRUE,rep(1,4))
check(qr(Hnear)$rank<2 && abs(out$se-.5)<1e-6 && is.finite(out$low) &&
 grepl('method=MASS::ginv',out$message,fixed=TRUE),'near-singular outcome bread uses ginv and records its method')
near$H<-matrix(1,2,2)
out<-maic_delta(near,rep(1,4),matrix(zz,4,1),0,matrix(.01),TRUE,rep(1,4))
check(abs(out$se-.5)<1e-10 && is.finite(out$low) &&
 grepl('ginv_residual_H=',out$message,fixed=TRUE),'exactly singular outcome bread permits an estimable target')
near$H<-matrix(1);near$score<-matrix(signs,4,1)
near$dp0<-matrix(0,4,1);near$dp1<-matrix(.2,4,1)
Znear<-cbind(zz,zz+1e-5*c(-1,1,-1,1))
out<-maic_delta(near,rep(1,4),Znear,c(0,0),diag(.01,2),TRUE,rep(1,4))
check(qr(crossprod(Znear)/4)$rank<2 && abs(out$se-.5)<1e-6 &&
 grepl('ginv_residual_Ag=',out$message,fixed=TRUE),'near-singular calibration bread returns SE and ginv diagnostics')
out<-maic_delta(near,rep(1,4),cbind(zz,zz),c(0,0),diag(.01,2),TRUE,rep(1,4))
check(abs(out$se-.5)<1e-10 && is.finite(out$low),
 'exactly singular calibration bread is handled by ginv')
probe<-suppressWarnings(run.anchored(joint_rr_gop_design(event='common',ess_ratio=.5),100,100,
 seed=2,score_spec=make_brm_score_spec(),run_exact=FALSE,
 fit_control=list(brm_max_step=1L,brm_optim_maxit=1L,fc_max_step=1L,
  brm_thres=1e-3,fc_thres=1e-3,brm_optim_reltol=2e-6)))
check(probe$fits$brm$optimization_input$thres==1e-3 &&
 probe$fits[['brm-FC']]$optimization_input$thres==1e-3,
 'custom tolerances reach the actual BRM and FC solver inputs')
check(probe$fits$brm$optimization_input$optim.reltol==2e-6 &&
 probe$fits$brm$optimization$optim_reltol==2e-6,
 'inner tolerance reaches the BRM diagnostics and saved retry inputs')
