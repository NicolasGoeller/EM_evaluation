
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(tictoc)
library(reshape2)
library(gridExtra)

#tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

theta_null = c(1)

k = length(theta0) 
q = length(theta_null)

#theta_start = theta0

# Load the log-likelihood, its derivative, and the hessian
source("assignment_2/Probit_LL.R")
source("assignment_2/Probit_LL_g.R")
source("assignment_2/Probit_LL_h.R")
source("assignment_2/Probit_LL_r.R")
source("assignment_2/Probit_J_1.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)
theta_hat_r_vec = matrix(0,num,k)

#J_1_inv = matrix(0,k,k)
#J_2_inv = matrix(0,k,k)

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
	result_r <- optim(par = theta_null, Probit_LL_r, y = y, x = x, theta0_r= 1, 
	                  method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
	theta_hat_r = c(1, result_r$par)
	
	theta_hat_r_vec[it,1:k] = theta_hat_r
	result_r$hessian
	## Wald
	C_theta = matrix(c(0,1),nrow=1,ncol=2)
	c_theta = (theta_hat_r[1]-theta_hat_r[2])
	V_n = solve(Probit_LL_h(y,x,theta_hat_r))#(1/n)*
	Wald = c_theta*solve(C_theta%*%V_n%*%t(C_theta))*c_theta
	#V_n = solve(C_theta%*%((1/n)*Probit_LL_h(y,x,theta_hat_r))%*%t(C_theta))
	#V_n = C_theta%*%(1/n)*(solve(Probit_LL_h(y,x,theta_hat_r)))%*%t(C_theta)

	#Wald = sum(theta_hat_r)*V_n*sum(theta_hat_r)
	
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

	# LR Probit_LL_r(y,x,theta_hat_r[2],1)
 	## Why is the bracket negative LL_r is larger than LL (even the same is LL(theta_r))
	#LR = 2*n*(Probit_LL(y,x,theta_hat)-Probit_LL_r(y,x,theta_hat_r[1],1))
 	LR = 2*n*(Probit_LL(y,x,theta_hat_r)-Probit_LL(y,x,theta_hat))
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

#test_data <- list("Wald"=Wald_vec, 
#                  "Score"=Score_vec, 
#                  "LR"=LR_vec) %>% 
#  as.data.frame() %>% 
#  pivot_longer(cols=1:3, names_to="test", values_to="test_stat")
#
#
#
#test_data %>% 
#  filter(test == "Score") %>% 
#  ggplot(aes(x=test_stat, fill=factor(test)))+
#  geom_histogram(aes(y = stat(count)/300),bins=30)+
#  #stat_function(fun=dchisq(x, df = k))#+
#  facet_wrap("test")+
#  labs(title = "Test statistic distributions (N=300)",
#       x= "Test statistic value", y="Count")+
#  theme_light()+
#  theme(legend.position = "none")
###ggsave("assignment_2/MLE_test_trinity_d.png", width=16,height=9,units = "cm")


png("assignment_2/wald_hist.png")
hist(Wald_vec,30, freq=FALSE, 
     xlab = "Wald test statistic", ylab = "Density", main="Distribution of Wald test statistic") 
curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
dev.off()

png("assignment_2/score_hist.png")
hist(Score_vec,30, freq=FALSE,
     xlab = "Score test statistic", ylab = "Density", main="Distribution of Score test statistic") 
curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
dev.off()

png("assignment_2/lr_hist.png")
hist(LR_vec,30, freq=FALSE,
     xlab = "LR test statistic", ylab = "Density", main="Distribution of LR test statistic") 
curve(dchisq(x, df = q), from = 0, to = 15, add=TRUE)
dev.off()

