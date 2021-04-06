# WVIC: association and interaction test for genetic markers and interval-censored survival outcomes

This tutorial demonstrates the implmentation of proposed interaction and association test for interval-censored survival outcomes, whose underlying model is semiparametric transformation model. The main function is _WVIC_comprisk_main.R_, which involves basic functions in _Base_Function_Transform_frailty.R_ and method for obtaining Nonparametric Maximum Likelihood Estimation (NPMLE) for regression coefficients and cumulative baseline function in _unpencoxIC.default.R_ (Li et. al. (2019) https://doi.org/10.1177/0962280219856238) or R functions start with _semiparam_transform_model_.

# Weighted V test (WVIC) usage
To perform the test, use function _WVIC_comprisk_main.R_ in 'R_code' file with the following input:
* X:           genetic variable to be tested
* Z:           adjustment covariates (default is NULL for no adjustment covariates)
* bound_times: 2-d matrix of (L,R) bound times (R=Inf if right censored)
* trunct:      left truncation times (default is NULL for no left truncation)
* Gsim:        Genetic similarity matrix
* Fsim:        Background similarity matrix (default is NULL for no heterogeneity effect)
* r:           semi-param transformation model parameter (r>=0)
* lim, acc:    parameters used in Davies method (default is lim=5e5, acc=1e-5)
* covmat:      calculate covariance matrix (default is FALSE)

output values are:
* coef_est: regression coefficient estimates
* CovEst: covariance matrix estimate
* pval: test p-value
* ifault: p-value calculation error index (refer to 'ifault' in Davies method in _CompQuadForm_ R package)
* lcr: left-censor rate
* rcr: right-censor rate

# Sample data analysis
Analysis of sample data is included in 
