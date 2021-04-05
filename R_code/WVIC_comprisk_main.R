
suppressMessages(library(MASS))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(ALassoSurvIC))
suppressMessages(library(CompQuadForm))

source("./unpencoxIC.default.R")
source("./Base_Function_Transform_frailty.R")
source("./semiparam_transform_model_rp_adj.R")
source("./semiparam_transform_model_rp_noadj.R")
source("./semiparam_transform_model_r0_adj.R")
source("./semiparam_transform_model_r0_noadj.R")


WVIC_test = function(X, Z=NULL, bound_times, trunct=NULL, Gsim, Fsim=NULL, r, lim=5*10^5, acc=10^(-5), covmat=FALSE){
  #-------------------------------------------------------------------------------#
  # X:           genetic variable to be tested
  # Z:           adjustment covariates (NULL if no adjustment covariates)
  # bound_times: 2-d matrix of (L,R) bound times (col2=Inf if right censored)
  # trunct:      left truncation times
  # Gsim:        Genetic similarity
  # Fsim:        Background similarity (NULL if no heterogeneity)
  # r:           semi-param transformation model parameter (r>=0)
  # lim, acc:    parameters used in Davies method
  # covmat:      TRUE or FALSE (calculate covariance matrix or not)
  #-------------------------------------------------------------------------------#
  
  if(r>0 & !is.null(trunct)){stop('truncation times not supported for semi-parametric transformation model with r>0')}
  
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  
  L_end = bound_times[,1]
  R_end = bound_times[,2]
  
  # left-, right-censor rates
  lc_rate = sum(L_end==0)/nrow(bound_times)
  rc_rate = sum(R_end==Inf)/nrow(bound_times)
  
  # Consider adjustment covariates
  if(!is.null(Z)){
    Z_name = colnames(as.data.frame(Z))
    Z = as.matrix(Z)
    if(r == 0){
      if(covmat){
        fit0 <- suppressMessages(unpencoxIC.default(L_end, R_end, Z, normalize.X=FALSE, covmat=TRUE))
        covar_matrix = fit0$cov
      } else {
        fit0 <- suppressMessages(unpencoxIC.default(L_end, R_end, Z, normalize.X=FALSE, covmat=FALSE))
        covar_matrix = NULL
      }
      jumppts = c(0,baseline(fit0)$upper.set)
      cumhaz = c(0,baseline(fit0)$clambda)
    }
    if(r > 0){
      if(covmat){
        fit0 <- suppressMessages(semipar_trans_fit_rp_adj(Z, bound_time=cbind(L_end,R_end), perturb=5, n.l=3, r_alpha=r, maxiter=100, cov_method='ZL1'))
        covar_matrix = fit0$covest
      } else {
        fit0 <- suppressMessages(semipar_trans_fit_rp_adj(Z, bound_time=cbind(L_end,R_end), perturb=5, n.l=3, r_alpha=r, maxiter=100, cov_method=0))
        covar_matrix = NULL
      }
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    regcoef = fit0$b
  }
  # Do not consider adjustment covariates
  if(is.null(Z)){
    if(r == 0){
      fit0 <- suppressMessages(semipar_trans_fit_r0_noadj(bound_time=cbind(L_end,R_end), r_alpha=r, maxiter=100))
      covar_matrix = NULL
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    if(r > 0){
      fit0 <- suppressMessages(semipar_trans_fit_rp_noadj(bound_time=cbind(L_end,R_end), n.l=3, r_alpha=r, maxiter=100))
      covar_matrix = NULL
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    regcoef = NULL
  }
  
  #message(" Now: calculating the p-value")
  
  if(is.null(Fsim)){
    R = Gsim
  } else {
    R = Gsim*(1+Fsim)
  }
  
  if(is.null(Z)){
    I = diag(n)
    J = matrix(1/n,n,n)
    R.res = (I-J) %*% R %*% (I-J)
  } else {
    I = diag(n)
    Z.bar = cbind(rep(1,n),Z)
    hatmatrix = Z.bar %*% solve(t(Z.bar)%*%Z.bar) %*% t(Z.bar)
    R.res = (I-hatmatrix) %*% R %*% (I-hatmatrix)
  }
  
  if(is.null(Z)){
    if(is.null(trunct)){
      s = apply(cbind(L_end,R_end), 1, FUN=s_cal_noadj, jumppts=jumppts, cumhaz=cumhaz, r=r)
    } else {
      s = apply(cbind(trunct,L_end,R_end), 1, FUN=s_cal_noadj_trunc, jumppts=jumppts, cumhaz=cumhaz, r=r)
    }
  } else {
    if(is.null(trunct)){
      s = apply(cbind(L_end,R_end,Z), 1, FUN=s_cal_adj, jumppts=jumppts, cumhaz=cumhaz, regcoef=regcoef, r=r)
    } else {
      s = apply(cbind(trunct,L_end,R_end,Z), 1, FUN=s_cal_adj_trunc, jumppts=jumppts, cumhaz=cumhaz, regcoef=regcoef, r=r)
    }
  }
  eta = eigen(R.res, symmetric=TRUE, only.values=TRUE)$values
  
  lambda1 = as.numeric(t(s) %*% s)
  nV = c(t(s) %*% R.res %*% s / n)
  
  if(is.null(Z)){
    Davies = davies(nV, lambda1*eta/(n^2), lim=lim, acc=acc)
  } else {
    Davies = davies(nV, lambda1*eta/(n*(n-ncol(Z)-1)), lim=lim, acc=acc)
  }
  
  output <- list()
  output$coef_est = regcoef         #regression coef estimates
  output$CovEst = covar_matrix      #covariance matrix estimate for regression coef
  output$pval = Davies$Qq           #p-value
  output$ifault = Davies$ifault     #p-value calculation error
  output$lcr = lc_rate              #left-censor rate
  output$rcr = rc_rate              #right-censor rate
  return(output)
} #end of function#


