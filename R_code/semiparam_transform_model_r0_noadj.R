
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
DentalEM_r0_noadj = function(bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(bound_time)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
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
    
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = rep(1,n)       # length n array
    
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
    
    ####################################
    ###         Lambda Update
    ####################################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_r0 function#

###################### Main function ######################
semipar_trans_fit_r0_noadj = function(bound_time, r_alpha, maxiter){
  #-----------------------------------------------------------
  # bound_time: left and right bound time
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_r0_noadj(bound_time, order_bound_times, r_alpha, maxiter)
  
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Output results ######################
  output = list()
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
} #end of output function#
