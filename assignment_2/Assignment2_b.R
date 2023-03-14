
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

theta_start = rep(0,k)

# Load the log-likelihood, its derivative, and the hessian
source("assignment_2/Probit_LL_h.R")
source("assignment_2/Probit_LL.R")
source("assignment_2/Probit_J_1.R")

num = 100			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)

J_1_inv = matrix(0,k,k)
J_2_inv = matrix(0,k,k)

B = 399

J_inv_boot = matrix(0,k,k)

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat = probit$coefficients
  #probit <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  #theta_hat = probit$par
  
  #switch out glm for our optim function 
  
  theta_hat_vec[it,1:k] = theta_hat
  
  J_1_inv = J_1_inv + solve(Probit_J_1(y,x,theta_hat))
  
  J_2_inv = J_2_inv + solve(Probit_LL_h(y,x,theta_hat))
  
  theta_hat_boot = matrix(0,B,k)
  
  for (b in 1:B) {
    index_b <- sample(length(y),length(y),replace=TRUE)
    y_b <- y[index_b]
    x_b <- x[index_b,]
    dat_b = data.frame(x_b,y_b)
    probit <- glm(y_b~x_b-1, data=dat_b, family = binomial(link = "probit"))
    #probit <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
    theta_hat_boot[b,1:k] = probit$coefficients
  }
  
  J_inv_boot = J_inv_boot + var(theta_hat_boot)
}

# Average of theta_hat over Monte Carlo iterations
colMeans(theta_hat_vec)

# Variance of theta_hat over Monte Carlo iterations
var(theta_hat_vec)

# Average of variance estimate based on J_1
J_1_inv/num

# Average of variance estimate based on J_2
J_2_inv/num

# Average of variance estimate based on J_3
J_inv_boot/num 

toc()


############ RESULTS!!!!############

#> # Average of theta_hat over Monte Carlo iterations
#  > colMeans(theta_hat_vec)
#[1] 1.0037004 0.9977425
#> 
#  > # Variance of theta_hat over Monte Carlo iterations
#  > var(theta_hat_vec)
#[,1]        [,2]
#[1,] 0.013638099 0.006865704
#[2,] 0.006865704 0.014315621
#> 
#  > # Average of variance estimate based on J_1
#  > J_1_inv/num
#[,1]        [,2]
#[1,] 0.012352861 0.007735274
#[2,] 0.007735274 0.017101566
#> 
#  > # Average of variance estimate based on J_2
#  > J_2_inv/num
#[,1]        [,2]
#[1,] 0.011931341 0.007104595
#[2,] 0.007104595 0.015858498
#> 
#  > # Average of variance estimate based on J_3
#  > J_3_inv/num
#[,1]        [,2]
#[1,] 0.011920912 0.007070086
#[2,] 0.007070086 0.015763052
#> 
#  > # Average of variance estimate based on J_3
#  > J_inv_boot/num 
#[,1]        [,2]
#[1,] 0.012534687 0.007547447
#[2,] 0.007547447 0.016350754
#> 
#  > toc()
#186.89 sec elapsed
#
#
#