library(tidyverse)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

# Load the log-likelihood, its derivative, and the hessian
source("assignment_1/Probit_LL.R")
source("assignment_1/Probit_NLS.R")
source("assignment_1/Probit_GMM.R")
source("assignment_1/Probit_J_1.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_ML_vec = matrix(0,num,k)
theta_hat_NLS_vec = matrix(0,num,k)
theta_hat_GMM_vec = matrix(0,num,k)

inside_N_ML = rep(0,num)
inside_N_NLS = rep(0,num)
inside_N_GMM = rep(0,num)

J_1_sum = matrix(rep(0,4),2)
J_2_sum = matrix(rep(0,4),2)
Var_hat_NLS_sum = matrix(rep(0,4),2)
Var_hat_GMM_sum = matrix(rep(0,4),2)

epsilon = 0.1

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  #u = rnorm(n,0,1)							 	# error term
  u = (dchisq(x,1)-1)/sqrt(2)  
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  # ML
  result <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"),
                  control = list(reltol=1e-9), hessian=TRUE) 
  
  theta_hat_ML = result$par
  theta_hat_ML_vec[it,1:k] = theta_hat_ML
  
  if (sqrt(sum((theta_hat_ML - theta0)^2)) < epsilon) {
    inside_N_ML[it] = 1
  }
  
  J_1_sum = J_1_sum + (1/n)*solve(Probit_J_1(y=y,x=x,par=result$par))
  J_2_sum = J_2_sum + (1/n)*solve(result$hessian)
  
  # NLS
  result <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"),
                  control = list(reltol=1e-9), hessian=TRUE)
  
  theta_hat_NLS = result$par
  theta_hat_NLS_vec[it,1:k] = theta_hat_NLS
  
  if (sqrt(sum((theta_hat_NLS - theta0)^2)) < epsilon) {
    inside_N_NLS[it] = 1
  }
  
  sigma_hat = Probit_Sigma_NLS(y=y,x=x,par=result$par)
  H_hat = solve((1/n)*result$hessian)
  Var_hat_NLS_sum = Var_hat_NLS_sum + (1/n)*(H_hat %*% sigma_hat %*% H_hat)
  
  # GMM
  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"),
                  control = list(reltol=1e-9), hessian=TRUE)
  
  theta_hat_GMM = result$par
  theta_hat_GMM_vec[it,1:k] = theta_hat_GMM
  
  if (sqrt(sum((theta_hat_GMM - theta0)^2)) < epsilon) {
    inside_N_GMM[it] = 1
  }
  
  #Var_hat_GMM_sum = Var_hat_GMM_sum + (1/n)*Probit_Var_GMM(y=y,x=x,par=result$par)
  Var_hat_GMM_sum = Var_hat_GMM_sum + (1/n)*Probit_Var_GMM(y=y,x=x,par=result$par)
}