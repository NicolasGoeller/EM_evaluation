library(tidyverse)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 3 			# Set the sample size

#theta0 = c(1,1) 	# True parameter value

#k = length(theta0) 

#source("assignment_1/Probit_LL.R")

b_val = c(10,100,1000)			# Number of Monte Carlo iterations

epsilon = 0.1
data <- c(-0.1,0.2,0.7)

bootstrap_df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(bootstrap_df) <- c("B","value")

for (num in b_val){
  b_est <- rep(0,num)
  for (it in 1:num) {
    # Data generating process
    x <- sample(data, size = n, rep=T) #Sampling 

    t <- mean(x) # Sample mean
    
    b_est[it] <- t
  }
  print("complete")
  b_data <- data.frame("B"=rep(num, num), "value"=b_est)
  bootstrap_df <- bind_rows(bootstrap_df, b_data)
  #bootstrap_list <- append(bootstrap_list, list(b_est))
  #mean_list <- c(mean_list, mean(b_est))
}
bootstrap_df %>% 
  mutate(B = paste0("B=",as.character(B))) %>% 
  ggplot(aes(x=value, fill=factor(B)))+
    stat_ecdf(geom = "step")+
    stat_ecdf(geom = "line", color="red")+
    #geom_smooth(method = glm(family = binomial((link="probit"))))
    facet_wrap("B")+
    labs(title = "Empirical Conditional Distribution Function",
       x= "Estimator value", y="")
ggsave("assignment_2/bootstrap_estim_a.png", width=16,height=9,units = "cm")
