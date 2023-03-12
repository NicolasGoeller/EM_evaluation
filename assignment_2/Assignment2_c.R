library(tidyverse)
library(margins)
library(gmm)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

# Load the log-likelihood, its derivative, and the hessian
source("assignment_2/Probit_LL.R")
source("assignment_2/Probit_NLS.R")
source("assignment_2/Probit_GMM.R")
source("assignment_2/avg_probit_effect.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_ML_vec = matrix(0,num,k)
theta_hat_NLS_vec = matrix(0,num,k)
theta_hat_GMM_vec = matrix(0,num,k)

avg_margeff_ML = rep(0,num)
avg_margeff_NLS = rep(0,num)
avg_margeff_GMM = rep(0,num)

inside_N_ML = rep(0,num)
inside_N_NLS = rep(0,num)
inside_N_GMM = rep(0,num)

epsilon = 0.1

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  #u = rnorm(n,0,1)							 	# error term
  u = (rchisq(n,1)-1)/sqrt(2)  
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  # ML
  #result <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"),
  #                control = list(reltol=1e-9), hessian=TRUE) 
  #theta_hat_ML = result$par
  
  dat = data.frame(x,y)
  result <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat_ML = result$coefficients
  
  theta_hat_ML_vec[it,1:k] = theta_hat_ML
  avg_margeff_ML[it] = avg_probit_effect(theta_hat_ML, 1, x)
  
  if (sqrt(sum((theta_hat_ML - theta0)^2)) < epsilon) {
    inside_N_ML[it] = 1
  }
  
  # NLS
  result <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"),
                  control = list(reltol=1e-9), hessian=TRUE)
  
  theta_hat_NLS = result$par
  theta_hat_NLS_vec[it,1:k] = theta_hat_NLS
  avg_margeff_NLS[it] = avg_probit_effect(theta_hat_NLS, 1, x)
  
  if (sqrt(sum((theta_hat_NLS - theta0)^2)) < epsilon) {
    inside_N_NLS[it] = 1
  }
  
  # GMM
  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"),
                  control = list(reltol=1e-9), hessian=TRUE)
  
  theta_hat_GMM = result$par
  theta_hat_GMM_vec[it,1:k] = theta_hat_GMM
  avg_margeff_GMM[it] = avg_probit_effect(theta_hat_GMM, 1, x)
  
  if (sqrt(sum((theta_hat_GMM - theta0)^2)) < epsilon) {
    inside_N_GMM[it] = 1
  }
}

#### Generate "True" sample

n_true = 10000

x_true = cbind(matrix(1,n_true,1),matrix(rnorm((k-1)*n_true,0,1),ncol=k-1)) 	# regressors

#u = rnorm(n,0,1)							 	# error term
u_true = (rchisq(n_true,1)-1)/sqrt(2)  

y_star_true = x_true %*% theta0 + u_true					# latent "utility"

y_true = ceiling(y_star_true/(max(abs(y_star_true))+0.1))				# observed outcome

# ML
dat_true = data.frame(x_true,y_true)
result <- glm(y_true~x_true-1, data=dat_true, family = binomial(link = "probit"))
theta_hat_ML = result$coefficients
true_ML_eff = avg_probit_effect(theta_hat_ML, 1, x_true)

# NLS
result <- optim(par = theta0, Probit_NLS, y = y_true, x = x_true, method = c("BFGS"),
                control = list(reltol=1e-9), hessian=TRUE)
theta_hat_NLS = result$par
true_NLS_eff = avg_probit_effect(theta_hat_NLS, 1, x_true)

# GMM
result <- optim(par = theta0, Probit_GMM, y = y_true, x = x_true, method = c("BFGS"),
                control = list(reltol=1e-9), hessian=TRUE)
theta_hat_GMM = result$par
true_GMM_eff = avg_probit_effect(theta_hat_GMM, 1, x_true)


# Plot the histogram of the distribution of the estimator over Monte Carlo iterations

theta1_data <- list("ML"=avg_margeff_ML, 
                    "NLS"=avg_margeff_NLS, 
                    "GMM"=avg_margeff_GMM) %>% 
  as.data.frame() %>% 
  pivot_longer(cols=1:3, names_to="type", values_to="estim")

theta1_data %>% 
  ggplot(aes(x=estim, fill=factor(type)))+
  geom_histogram(bins=30)+
  geom_vline(data = data.frame(xint=true_ML_eff,type="ML"), aes(xintercept = xint))+#, linetype = "line")
  geom_vline(data = data.frame(xint=true_NLS_eff,type="NLS"), aes(xintercept = xint))+#, linetype = "line")
  geom_vline(data = data.frame(xint=true_GMM_eff,type="GMM"), aes(xintercept = xint))+#, linetype = "line")
  facet_wrap("type")+
  labs(title = "Avg. Marginal Effects for Different Estimators (N=300)",
       subtitle = "Black vertical lines denote 'true' value from sample N=10000",
       x= "Avg Marginal Effect", y="Count")+
  theme_light()+
  theme(legend.position = "none")
#ggsave("assignment_2/avg_marg_eff_estim_c.png", width=16,height=9,units = "cm")
