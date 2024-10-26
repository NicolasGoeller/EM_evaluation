library(tidyverse)
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

## Load function for MC Probit
source("assignment_final/probit_beta_test.R")

## Set DGP parameters
beta = seq(-2,2,0.1)	### Build in constraints for estimated ranges	
sigmasq = c(0,0.05,0.15)#^2
dgp_param <- merge(beta,sigmasq) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.list()

## Set MC parameters
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
  print(param)
  
  ## Run individual MC iterations
  test_data = probit_beta_test(n, num,param, theta_null)
  final_data = rbind(final_data,test_data)

} 

write_csv(final_data, "assignment_final/beta_test_data.csv")
toc()

plot <- final_data %>% 
  rename(Wald = wald,Score = score,`Likelihood Ratio`=lr) %>% 
  pivot_longer(cols= c(Wald,Score,`Likelihood Ratio`),names_to = "test_stat",values_to = "test_value") %>% 
  mutate(reject_test = if_else(test_value > cv,1,0),
         sigma_star = paste0("Sigma*: ",sigma_star)) %>% 
  group_by(beta_star,sigma_star,test_stat) %>% 
  summarise(mean_rejection = mean(reject_test,na.rm = T)) %>% 
  ggplot(aes(x=beta_star,y=mean_rejection, color=test_stat))+
  geom_line()+
  facet_wrap(~sigma_star)+
  labs(x="True Mean of Slope distribution (Beta*)",y="Mean Rejection Probability",
       color='Test Statistic')
ggsave(plot,"assignment_final/beta_test_v1-1.png",width = 16, height=9, units="cm")


