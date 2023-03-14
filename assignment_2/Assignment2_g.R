
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(tictoc)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

theta_null = c(1,1)

k = length(theta0) 
q = length(theta_null)

# Load the log-likelihood, its derivative, and the hessian
source("assignment_2/Probit_GMM.R")
source("assignment_2/Probit_GMM_r.R")
source("assignment_2/Probit_Var_GMM.R")
source("assignment_2/GMM_Omega.R")
source("assignment_2/GMM_QLR.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)
theta_hat_r_vec = matrix(0,num,k)

Wald_vec = rep(0,num)
QLR_vec = rep(0,num)

reject_Wald = rep(0,num)
reject_Score = rep(0,num)
reject_QLR = rep(0,num)

cv = qchisq(.95, df=q) 
print("start")
for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  ### Constrained model optimisation
  result_r <- optim(par = theta_null[1], Probit_GMM_r, y = y, x = x, theta1_r= theta_null[2], 
                    method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_r = c(result_r$par, theta_null[2])
  
  theta_hat_r_vec[it,1:k] = theta_hat_r

  ## Wald
  R = matrix(c(1,-1),nrow=1,ncol=2)
  con = R%*%theta_hat_r - R%*%theta_null

  V_n = Probit_Var_GMM(y,x,theta_hat_r)
  Wald = con*solve(R%*%V_n%*%t(R))*con
  Wald_vec[it] = Wald
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1
    }
  }else{
    print("NA")
  }
  
  ## Unconstrained model optimisation
  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), 
                  control = list(reltol=1e-9), hessian=TRUE)
  theta_hat = result$par
  
  theta_hat_vec[it,1:k] = theta_hat
  
  # LR
  QLR = GMM_QLR(y,x,theta_hat_r) - GMM_QLR(y,x,theta_hat)
  QLR_vec[it] = QLR
  
  if (is.nan(QLR) == 0) {
    if (QLR > cv) {
      reject_QLR[it] = 1
    }
  }else{
    print("NA")
  }
}

mean(reject_Wald)
mean(reject_QLR)
