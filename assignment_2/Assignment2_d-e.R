
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(reshape2)
library(gridExtra)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

theta_null = c(1,1)
#theta_null = c(1,1.1)       
#theta_null = c(1,0.9)

k = length(theta0) 
q = 1

# Load the log-likelihood, its derivative, and the hessian
source("assignment_2/Probit_LL.R")
source("assignment_2/Probit_LL_g.R")
source("assignment_2/Probit_LL_h.R")
source("assignment_2/Probit_LL_r.R")
source("assignment_2/Probit_J_1.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)
theta_hat_r_vec = matrix(0,num,k)

Wald_vec = rep(0,num)
Score_vec = rep(0,num)
LR_vec = rep(0,num)

reject_Wald = rep(0,num)
reject_Score = rep(0,num)
reject_LR = rep(0,num)

cv = qchisq(.95, df=q) 

for (it in 1:num) {
	 # Data generating process
	 x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
 
	 u = rnorm(n,0,1)							 	# error term
 
	 y_star = x %*% theta0 + u						# latent "utility"
 
	 y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
 
	 ### Constrained model optimisation
	 result_r <- optim(par = theta_null[1], Probit_LL_r, y = y, x = x, theta1_r= theta_null[2], 
                     method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
   theta_hat_r = c(result_r$par, theta_null[2])
   
   theta_hat_r_vec[it,1:k] = theta_hat_r
 
   ## Wald
   R = matrix(c(1,-1),nrow=1,ncol=2)
   con = R %*% theta_hat_r - R%*%theta_null
   V_n = solve(Probit_LL_h(y,x,theta_hat_r))
   
   Wald = con*solve(R%*%V_n%*%t(R))*con
   
   Wald_vec[it] = Wald
   
   if (is.nan(Wald) == 0) {
   	if (Wald > cv) {
   	  reject_Wald[it] = 1
   	}
   }else{
     print("NA")
   }
   
   ## Score
   Score = Probit_LL_g(y,x,theta_hat_r)
  
   Score_stat = t(Score)%*%solve((1/n)*Probit_J_1(y,x,theta_hat_r))%*%Score
   Score_vec[it] = Score_stat
   
   if (is.nan(Score_stat) == 0) {
   	if (Score_stat > cv) {
   	  reject_Score[it] = 1
   	}
   }else{
     print("NA")
   }
   
   ## Unconstrained model optimisation
   result <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), 
                   control = list(reltol=1e-9), hessian=TRUE)
   theta_hat = result$par
   
   theta_hat_vec[it,1:k] = theta_hat

   # LR
   LR = 2*n*(Probit_LL(y,x,theta_hat_r)-Probit_LL(y,x,theta_hat_r))
   LR_vec[it] = LR
   
   if (is.nan(LR) == 0) {
     if (LR > cv) {
       reject_LR[it] = 1
     }
   }else{
     print("NA")
   }
} 

mean(reject_Wald)
mean(reject_Score)
mean(reject_LR)

## test theta1 = 1
#> mean(reject_Wald)
#> mean(reject_Wald)
#[1] 0.015
#> mean(reject_Score)
#[1] 0.065
#> mean(reject_LR)
#[1] 0.045

## test theta1 = 1.1
#> mean(reject_Wald)
#[1] 0.038
#> mean(reject_Score)
#[1] 0.072
#> mean(reject_LR)
#[1] 0.121

## test theta1 = 0.9
#> mean(reject_Wald)
#[1] 0.018
#> mean(reject_Score)
#[1] 0.221
#> mean(reject_LR)
#[1] 0.151

#### Histograms for test statistics
#png("assignment_2/wald_hist.png")
#hist(Wald_vec,30, freq=FALSE, 
#     xlab = "Wald test statistic", ylab = "Density", main="Distribution of Wald test statistic") 
#curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
#dev.off()
#
#png("assignment_2/score_hist.png")
#hist(Score_vec,30, freq=FALSE,
#     xlab = "Score test statistic", ylab = "Density", main="Distribution of Score test statistic") 
#curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
#dev.off()
#
#png("assignment_2/lr_hist.png")
#hist(LR_vec,30, freq=FALSE,
#     xlab = "LR test statistic", ylab = "Density", main="Distribution of LR test statistic") 
#curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
#dev.off()
