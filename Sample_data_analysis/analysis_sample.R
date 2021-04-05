
load('sample_data.RData')

setwd('../R_code/')
source('WVIC_comprisk_main.R')

test_results = WVIC_test(X=X, Z=Z, bound_times=LR, trunct=trunct, Gsim=Gsim, Fsim=NULL, r=r, lim=5*10^5, acc=10^(-5), covmat=FALSE)
pval = test_results$pval
