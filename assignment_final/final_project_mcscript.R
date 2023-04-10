library(tidyverse)

rm(list=ls()) 		# Clear workspace

## Set DGP parameters
beta1 = seq(0,2,0.1)			
sigma = c(0,0.05,0.15)
dgp_param <- merge(beta1,sigma) %>% t() %>% as.data.frame() %>% as.list()

## Set MC parameters
epsilon = 0.1     # Set error margin
n = 300 			    # Set the sample size
num = 1000        # Set number of MC runs

## Create final data set
columns <- c("beta1","sigma","wald","score","lr")
final_data <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(final_data) <- columns


for (param in dgp_param){
  

  #Random coefficient model: beta0, beta1, sigma
  theta0 = list(mean= c(1,param[1]),
                var= c(0,param[2])) 
  
  test_data <- mcprobit_beta1_test(n, num, param)
  
  final_data <- rbind(final_data,test_data)
  
}





####First case: plot test values (y-axis) and mean beta1 values (x-axis) over the three tests and with different beta1 variances
#### Implies #MCdraws(num)*#va-values(k)*#beta1-options(i) - each column of 3xk graph is one MC run for a var value

####Second case: plot test values (y-axis) and 