
suppressMessages(library(MASS))
suppressMessages(library(statmod))
suppressMessages(library(nlme))
#suppressMessages(library(survival))
#suppressMessages(library(Rcpp))
#suppressMessages(library(RcppArmadillo))

Gfunc = function(x, r){
  if(r==0){ x }
  else if(r>0){ log(1+r*x)/r }
}

###################### Function for semi-param transformation model fit ######################
DentalEM_r0_adj = function(X, bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(X)
  p = ncol(X)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  gamma.ini = rep(0, p)
  
  Params_old = rep(0, p)
  Params_new = gamma.ini
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  #-------------------- Apply EM Algorithm --------------------#
  ### Initial values
  iter = 0
  
  ### collect all loglik values
  loglik_all = NULL
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  ### Start iteration
  while (iter < maxiter){
    
    iter = iter + 1
    
    Params_old = Params_new
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    gamma = as.matrix(Params_old)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = c(exp(X %*% gamma))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1)    #n by 1 matrix
    
    ### logLikelihood
    lik_i = rep(NA,n)
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    loglik_i = log(lik_i)
    loglik = sum(loglik_i)
    loglik_all = c(loglik_all, loglik)
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = ((1/(1-exp(-sum.lambda.exp.between[i,])))*exp.term[i]*lambda_old*jumptidx_in_LR[,i])
      }
    }
    
    ### Score
    sum.exp.term = c(outer(order_bound_times, Rij.star, '<=') %*% matrix(exp.term, ncol=1))
    
    Score_vec = NULL
    for(j in 1:p){
      score.step1 = sapply(1:length(order_bound_times), function(y){ sum((order_bound_times[y]<=Rij.star) * X[,j] * exp.term) })
      score = sum(sapply(1:n, function(x){
        sum((order_bound_times<=Rij.star[x]) * X[x,j] * E.Cijt[x,]) - sum((order_bound_times<=Rij.star[x]) * E.Cijt[x,] * score.step1/sum.exp.term)
      }))
      Score_vec = rbind(Score_vec, score)
    }
    
    ### Hessian
    sum.X.exp.term = NULL
    for(j in 1:p){
      sum.X.exp.term = rbind(sum.X.exp.term, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*exp.term) }))
    }
    
    sum.XX.exp.term = NULL
    for(j in 1:p){
      for(k in j:p){
        sum.XX.exp.term = rbind(sum.XX.exp.term, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*X[,k]*exp.term) }))
      }
    }
    
    Hessian = matrix(0, nrow=p, ncol=p)
    # X single columns
    idx = 0
    for(i in 1:p){
      for(j in i:p){
        idx = idx + 1
        Hessian[i,j] = sum(sapply(1:n, function(x){
          sum(E.Cijt[x,]*((sum.XX.exp.term[idx,]*sum.exp.term - sum.X.exp.term[i,]*sum.X.exp.term[j,]) / (sum.exp.term)^2)*(order_bound_times<=Rij.star[x]))
        }))
      }
    }
    
    #Hessian_mat = Hessian
    Hessian_mat = Hessian + t(Hessian)
    diag(Hessian_mat) = 0.5*diag(Hessian_mat)
    
    ### Update parameters
    #Params_new = Params_old + c(solve(Hessian_mat) %*% Score_vec)
    Params_new = Params_old + c(ginv(Hessian_mat) %*% Score_vec)
    
    ####################################
    ###         Lambda Update
    ####################################
    gamma.l = as.matrix(Params_new)
    exp.term.l = c(exp(X %*% gamma.l))       # length n array
    
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term.l * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.param = max(abs(c(Params_new - Params_old)))
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(c(error.param, error.lambda))<=0.0001){ break }
    #if(max(error.param) <= 0.0005){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("Param change: %0.6f",error.param))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    # newpar = as.data.frame(matrix(Params_new, nrow=1))
    # colnames(newpar) = c(paste0('gamma',c(1:p)))
    # print(newpar)
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['param']] = Params_new
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_r0_adj function#

###################### Function for semi-param transformation model fit (w/ fixed gamma) ######################
dentalEM_fixgamma_r0_adj = function(X, gamma_est, bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  # gamma_est: gamma estimates
  #-----------------------------------------------------------
  n = nrow(X)
  p = ncol(X)
  ###################### Initial Values ######################
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  ###################### Apply EM Algorithm ######################
  iter = 0
  while (iter < maxiter){
    
    iter = iter + 1
    lambda_old = lambda_new
    gamma = as.matrix(gamma_est)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = c(exp(X %*% gamma))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1)     #n by n.l matrix
    
    ### logLikelihood
    lik_i = rep(NA,n)
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    loglik_i = log(lik_i)
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = ((1/(1-exp(-sum.lambda.exp.between[i,])))*exp.term[i]*lambda_old*jumptidx_in_LR[,i])
      }
    }
    
    ###################### Lambda Update ######################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ###################### Stopping Criteria ######################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
  } # end of while loop #
  
  return(loglik_i)
  
} #end of DentalEM_r0_adj (fixed gamma) function#

###################### Main function ######################
semipar_trans_fit_r0_adj = function(X, bound_time, perturb, r_alpha, maxiter, cov_method){
  #-----------------------------------------------------------
  # X: covariate matrix
  # bound_time: left and right bound time
  # perturb: integer 1, 5, 10, for perturbation
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  # cov_method = 0: do not estimate covariance_mat
  #            = integer: boostrap method for covariance_mat
  #            = 'ZL1': 1st-order (Zeng 2017)
  #            = 'ZL2': 2nd-order (Zeng 2016)
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  n = nrow(X)
  p = ncol(X)
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_r0_adj(X, bound_time, order_bound_times, r_alpha, maxiter)
  
  param_bl = fit_bl$param          #parameter estimate
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Zeng 2017 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL1'){
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    param_delta = matrix(rep(param_bl,p),byrow=T,nrow=p) + diag(delta,p,p)
    
    #Gradient of logliklihood
    loglik_grad = NULL
    loglik_0 = dentalEM_fixgamma_r0_adj(X, param_bl, bound_time, order_bound_times, r_alpha, maxiter)
    for(i in 1:nrow(param_delta)){
      loglik_delta = dentalEM_fixgamma_r0_adj(X, param_delta[i,], bound_time, order_bound_times, r_alpha, maxiter)
      loglik_grad = rbind(loglik_grad, (loglik_delta-loglik_0)/delta)
    }
    
    # V matrix
    V_mat = 0
    for(j in 1:n){ V_mat = V_mat + tcrossprod(as.matrix(loglik_grad[,j])) }
    
    # covariance matrix estimate
    covest = solve(V_mat)
    
  } #end of if(cov_method=='ZL1')#
  ###################### Zeng 2016 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL2'){
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    
    #Gradient of logliklihood
    V_mat = matrix(NA, nrow=p, nrow=p)
    for(i in 1:p){
      for(j in i:p){
        param_delta_i = param_bl + delta*(c(1:p)==i)
        param_delta_j = param_bl + delta*(c(1:p)==j)
        param_delta_ij = param_bl + delta*(c(1:p)==i) + delta*(c(1:p)==j)
        V_mat[i,j] = sum(dentalEM_fixgamma_r0_adj(X, param_bl, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_i, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_j, bound_time, order_bound_times, n.l, r_alpha, maxiter)) +
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_ij, bound_time, order_bound_times, n.l, r_alpha, maxiter))
      }
    }
    V_mat = V_mat + t(V_mat)
    diag(V_mat) = 0.5*diag(V_mat)
    
    # covariance matrix estimate
    covest = -1*solve(V_mat)
    
  } #end of if(cov_method=='ZL2')#
  ###################### Do not need covariance matrix estimate ######################
  if(cov_method == 0){ covest = NULL }
  
  ###################### Output results ######################
  output = list()
  output[['b']] = c(param_bl)
  output[['covest']] = covest
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
  
} #end of output function#
