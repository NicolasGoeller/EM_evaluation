library(tidyverse)
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

## Load function for MC Probit
source("assignment_final/probit_beta1_test.R")

## Set DGP parameters
beta = seq(-2,2,0.1)	### Build in constraints for estimated ranges	
sigmasq = c(0,0.05,0.15)#^2
dgp_param <- merge(beta,sigmasq) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.list()

## Set MC parameters
epsilon = 0.1     # Set error margin
n = 300 			    # Set the sample size
num = 1000       # Set number of MC runs
cv = qchisq(.95, df=1) 

## Create final data set
columns <- c("beta_star","sigma_star","wald","score","lr","beta",  
             "sigma","beta_c","sigma_c")
final_data <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(final_data) <- columns

## Set Hypothesis values
theta_null <- c(0) #beta


for (param in dgp_param){
  set.seed(1) #Seed for random number generator
  i = 1
  print(param)
  
  while (i < num){
  
    if (i %in% c(1,500,1000)){
      print(i)
    }
    ## Run individual MC iterations and return NULL for error
    test_data <- tryCatch({
      probit_beta1_test(n, param, theta_null)
      },
      error= function(e){
        return(NULL)
        })

    ## If result is NULL rerun iteration with next draw
    if (is.null(test_data) == TRUE){
      print("Repeat")
      next
    ## Otherwise add to table and move to next parameter
    }else{ #
      final_data <- rbind(final_data,test_data)
      i = i +1
    }
  } 
}
write_csv(final_data, "assignment_final/beta_test_data.csv")
toc()

final_data %>% 
  pivot_longer(cols= c(wald,score,lr),names_to = "test_stat",values_to = "test_value") %>% 
  mutate(reject_test = if_else(test_value > cv,1,0)) %>% 
  group_by(beta_star,sigma_star,test_stat) %>% 
  summarise(mean_rejection = mean(reject_test,na.rm = T)) %>% 
  ggplot(aes(x=beta_star,y=mean_rejection, color=test_stat))+
  geom_line()+
  facet_grid(test_stat~sigma_star)

####First case: plot test values (y-axis) and mean beta1 values (x-axis) over the three tests and with different beta1 variances
#### Implies #MCdraws(num)*#va-values(k)*#beta1-options(i) - each column of 3xk graph is one MC run for a var value

####Second case: plot test values (y-axis) and 