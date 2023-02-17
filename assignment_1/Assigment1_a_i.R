
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

# Data generating process
x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors

u = rnorm(n,0,1)							 	# error term

y_star = x %*% theta0 + u						# latent "utility"

y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome

# Load the log-likelihood
source("assignment_1/Probit_LL.R")

# optim without user-specified gradient
result_b <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_b$par

# Load the derivative of the log-likelihood
source("assignment_1/Probit_LL_g.R")

# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_LL_g(y,x,theta0)
grad(function(u) Probit_LL(y,x,u),theta0)

# optim with user-specified gradient
result_c <- optim(par = theta0, Probit_LL, y = y, x = x, gr = Probit_LL_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_c$par

# Load the objective function for NLS
source("assignment_1/Probit_NLS.R")

# optim without user-specified gradient
result_e <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_e$par

# Load the derivative of the NLS objective function
source("assignment_1/Probit_NLS_g.R")

# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_NLS_g(y,x,theta0)
grad(function(u) Probit_NLS(y,x,u),theta0)

# optim with user-specified gradient
result_f <- optim(par = theta0, Probit_NLS, y = y, x = x, gr = Probit_NLS_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_f$par

# Load the objective function for MM
source("assignment_1/Probit_GMM.R")

# optim without user-specified gradient
result_h <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_h$par

# Load the derivative of the GMM objective function
source("assignment_1/Probit_GMM_g.R")
# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_GMM_g(y,x,theta0)
grad(function(u) Probit_GMM(y,x,u),theta0)

# optim with user-specified gradient
result_i <- optim(par = theta0, Probit_GMM, y = y, x = x, gr = Probit_GMM_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_i$par




