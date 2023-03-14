
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1,1) 	# True parameter value

theta_null = c(1,1,1) 

k = length(theta0) 

theta_start = rep(0,k)

# Load the log-likelihood, its derivative, and the hessian
source("Probit_LL.R")
source("Probit_LL_g.R")
source("Probit_LL_h.R")
source("J_1.R")
source("J_3.R")

num = 1000			# Number of Monte Carlo iterations

reject_Score = matrix(0,num,1)
reject_Score_boot = matrix(0,num,1)

B = 399

cv = qchisq(.95, df=3) 

for (it in 1:num) {
	# Data generating process
	x = cbind(matrix(1,n,1),matrix(rnorm(2*n,0,1),ncol=2)) 	# regressors

	u = rnorm(n,0,1)							 	# error term

	y_star = x %*% theta0 + u						# latent "utility"

	y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
	
	dat = data.frame(x,y)
	probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
	theta_hat = probit$coefficients
	
	# Score
	Score = Probit_LL_g(y,x,theta_null)

	Score = t(Score)%*%solve(J_1(y,x,theta_null))%*%Score
	if (is.nan(Score) == 0) {
  	if (Score > cv) {
  	  reject_Score[it] = 1
  	}
	}
	
	Score_boot = rep(0,B)
	
	for (b in 1:B) {
	  index_b <- sample(length(y),length(y),replace=TRUE)
	  y_b <- y[index_b]
	  x_b <- x[index_b,]
	  Score_b = Probit_LL_g(y_b,x_b,theta_hat)
	  Score_b = t(Score_b)%*%solve(J_1(y_b,x_b,theta_hat))%*%Score_b
	  Score_boot[b] = Score_b
	}

	Score_boot <- sort(Score_boot)
	cv_b <- Score_boot[floor(B*0.95)]

	if (is.nan(Score) == 0) {
	  if (Score > cv_b) {
	    reject_Score_boot[it] = 1
	  }
	}
}

mean(reject_Score)
mean(reject_Score_boot)

toc()


