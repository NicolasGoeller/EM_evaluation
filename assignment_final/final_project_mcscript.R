library(tidyverse)
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

## Load function for MC Probit
source("assignment_final/mcprobit_beta1_test.R")

## Set DGP parameters
beta0 = c(1)
beta1 = seq(-2,2,0.1)	### maybe constrain more at 2,-2 doesnt run, does for 1.5		
sigma = c(0,0.05,0.15)
dgp_param <- merge(beta1,sigma) %>% 
  cbind("z"=rep(1,length(beta1)*length(sigma)),.) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.list()

## Set MC parameters
epsilon = 0.1     # Set error margin
n = 300 			    # Set the sample size
num = 1000        # Set number of MC runs

## Create final data set
columns <- c("beta1","sigma","wald","score","lr")
final_data <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(final_data) <- columns

## Set Hypothesis values
theta_null <- c(1,0,0.05) #beta0,beta1,sigma


for (param in dgp_param){
  print(param)
  
  test_data <- mcprobit_beta1_test(n, num, param, theta_null)
  
  final_data <- rbind(final_data,test_data)
  
}

toc()



####First case: plot test values (y-axis) and mean beta1 values (x-axis) over the three tests and with different beta1 variances
#### Implies #MCdraws(num)*#va-values(k)*#beta1-options(i) - each column of 3xk graph is one MC run for a var value

####Second case: plot test values (y-axis) and 