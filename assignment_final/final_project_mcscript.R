library(tidyverse)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = list(mean= c(1,1),
              var= c(0,1)) #Random coefficient model: beta0, beta1, sigma

b_val = seq(0,2,0.1)			# Number of Monte Carlo iterations

epsilon = 0.1

for (num in b_val){
  set.seed(1)

  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1 + (x**2)%*%(theta0[[2]]**2))				# error term
  
  y_star = x %*% theta0[[1]] + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))
}