# This is a tutorial for the simulation process for Paper "A GENERAL THEORY OF HYPOTHESIS TESTS AND CONFIDENCE REGIONS FOR SPARSE HIGH DIMENSIONAL MODELS". 
# Model: Linear Regression with Decorrelated Score Method 
# Main settings: d (dimension), rho, n = 500;
# Throughout the simulation study, I first set the data generator process (DGP) of the covariates X:
# n = 100 independent and identical distribution samples with a multivariate Gaussian distribution
# Nd(0,Σ), where d = 100,200,500 and Σ is a diagonal-constant matrix with Σij = ρ|i−j|. ρ has
# four potential values, namely, 0.25, 0.4, 0.6, and 0.75. 

library(MASS)
# You can choose "glmnet" for lasso or "ncvreg" for SCAD.
# You can choose a family type to make lienar, logistic, and Poisson regression
library(glmnet)
library(ncvreg)
library(hdme)
Simulation<-function(n,d,rho){
  Sigma <- array(1,dim = c(d, d))
  for (i in 1:d){
    for (j in 1:d){
      Sigma[i,j]=rho^(abs(i-j))
    }
  }
  
  mean_<- array(0, dim = c(d,1))
  x = mvrnorm(n,mean_,Sigma) # DGP
  # Generate Y and beta
  beta <- array(0, dim = c(d,1))
  beta[2] = 1
  beta[3] = 1
  beta[4] = 1
  # You can use this code to implement linear regression, Logistic regress, and Poisson regression.
  #Linear regression
  #y = x %*% beta
  #Poisson regression
  # y <- rpois(n, exp(x %*% beta))
  #Logistic regression, see below.
  y <- rbinom(n, 1, (1 + exp(-x %*% beta))^(-1))
  fitcv<-cv.glmnet(x, y, family="binomial", alpha=1)
  beta_est = coef(fitcv, fitcv$lambda.min)
  # "Glmnet" is a package that fits generalized linear and similar models via penalized maximum likelihood.
  # You can choose a family type to make linear, logistic, and Poisson regression
  # For example, fitcv<-cv.ncvreg(x, y, penalty='SCAD',family="poisson" ).
  # fitcv<-ncvreg(x, y, penalty='SCAD',family="poisson")
  # beta_est = coef(fitcv, fitcv$lambda.min)
  
  c = array(0, dim = c(500, 1))
  # This part is to estimate a Dantzig Selector. 
  # For linear regression, we do not need to calculate x_,z_
  # x_,z_ are to calculate the Dantzig Selector. 
  # For Poisson regression, x_,z_ is set based on partial Fisher information. See paper Ning, 2017.
  for (i in 1:n ) {
    c[i,1] =  sqrt(exp( x[i, 1:d]%*%beta_est[2:(d+1) ] ) / (1 + exp( x[i, 1:d]%*%beta_est[2:(d+1) ] ) )^2 )
  }
  c_ = sqrt(c)
  x_ = x
  for (i in n){
    x_[i, 1:d] = x_[i,1:d]*c[i,1]
  }
  
  z_ = x[1:n,1]
  x_ = x[1:n, 2:d]
  # Cross-validation
  # Generalized Dantzig Selector with cross-validation.
  cv_fit <- cv_gds(x_, z_, family = "gaussian", no_lambda = 50, n_folds = 10)
  fit_gds = gds(x_, z_, family = "gaussian", lambda = cv_fit$lambda_min)
  gamma = coef(fit_gds)
  
  #calculate sum of score
  score=0
  score_005 = 0
  score_010 = 0
  score_015 = 0
  score_020 = 0
  score_025 = 0
  score_030 = 0
  score_035 = 0
  score_040 = 0
  score_045 = 0
  score_050 = 0
  score_055 = 0
  score_interval = 0
  I=0 # Partial Fisher Information
  for (j_2 in (1:n)){
    wx=0
    for(j_3 in 1:dim(gamma)[1]){
      wx = wx + x_[j_2, gamma[j_3,1]]*gamma[(j_3),2]
    }
    if(dim(gamma)[1] == 0  ){
      wx = 0
    }
    
    score_interval = score_interval + (y[j_2] - (exp(x[j_2, 1:d]%*%beta_est[2:(d+1) ])
                                                 /(1+exp(beta_est[1]+ x[j_2, 1:d]%*%beta_est[2:(d+1)])) )   )*(x[j_2, 1]-wx)
    score=score+ (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)])
                            /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)])) )   )*(x[j_2, 1]-wx)
    score_005 = score_005 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.05 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.05*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_010 = score_010 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.10 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.10*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_015 = score_015 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.15 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.15*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_020 = score_020 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.20 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.20*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_025 = score_025 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.25 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.25*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_030 = score_030 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.30 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.30*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_035 = score_035 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.35 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.35*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_040 = score_040 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.40 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.40*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_045 = score_045 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.45 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.45*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_050 = score_050 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.50 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.50*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    score_055 = score_055 + (y[j_2] - (exp(x[j_2, 2:d]%*%beta_est[3:(d+1)]  + 0.55 * x[j_2, 1])
                                       /(1+exp( x[j_2, 2:d]%*%beta_est[3:(d+1)] + 0.55*x[j_2, 1]  ) ) )   )*(x[j_2, 1]-wx)
    # Partial Fisher Information
    I = I + (exp(x[j_2, 1:d]%*%beta_est[2:(d+1)]) 
             / (1+exp(x[j_2, 1:d]%*%beta_est[2:(d+1)]) )^2) *x[j_2,1]*(x[j_2,1] - wx) 
  }  
  
  I = abs(I)/n
  score = (-1/n)*score
  U = sqrt(n)*score*I^(-0.5)
  return(U)
}


# Main Function
main<-function(){
  # Set n, d, rho;
  n = 500
  d = 100
  rho = 0.25
  result <- array(0, dim = c(500,1))
  T = 500 # Simulation times = 500;
  for (m in 1:T){
    result[m,1] = Simulation(n, d, rho)
  }
}

main()
