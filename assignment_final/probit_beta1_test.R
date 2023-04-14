probit_beta1_test <- function(n,param,theta_null){
  
  #set.seed(1) #Seed for random number generator
  
  ## Load functions
  source("assignment_final/rnd_Probit_LL.R")
  source("assignment_final/rnd_Probit_LL_cons.R")
  source("assignment_final/rnd_Probit_LL_g.R")
  source("assignment_final/rnd_Probit_J_1.R")
 
  line = 1
  test_data <- data.frame("beta_star"=rep(param[1],line),"sigma_star"=rep(param[2],line),
                          "wald"=rep(NA,line),"score"=rep(NA,line),"lr"=rep(NA,line),
                          "beta"=rep(NA,line),"sigma"=rep(NA,line),
                          "beta_c"=rep(NA,line),"sigma_c"=rep(NA,line))

  k = length(param)
  J_1_inv_sum = matrix(c(0,0,0,0),ncol=k)
  J_2_inv_sum = matrix(c(0,0,0,0),ncol=k)
  
#  for (i in 1:num){
    #print(i)
    ## Random Coefficient Probit DGP
    x = matrix(rnorm((k-1)*n,0,1),ncol=k-1) 	# regressors

    u = rnorm(n,0,1)                            # error term 1
    e = rnorm(n,0,1)                            # error term 2
    beta_i = param[1]+param[2]*e                # varying coefficients
    y_star = x*beta_i + u                       # latent "utility"

    #u = rnorm(n,mean=0, sd=sqrt(1+(t(x)%*%x)*param[2]^2))				# error term
    #y_star = x * param[1] + u						# latent "utility"
    y = ceiling(y_star/(max(abs(y_star))+0.1))

    ## Estimate unrestricted model
    result <- optim(par = param, rnd_Probit_LL, y = y, x = x,
                    method = c("L-BFGS-B"), lower= c(-Inf,0),
                    control = list(pgtol=1e-9), hessian=TRUE)

    theta_hat_ML = result$par
    test_data[line,"beta"]= theta_hat_ML[1]
    test_data[line,"sigma"]= theta_hat_ML[2]
    

    ## Wald Test
    J_2_inv_hat = solve(result$hessian)/n
    wald = (theta_hat_ML[1] - theta_null)^2/J_2_inv_hat[1,1]
    
    test_data[line,"wald"]= wald
    J_2_inv_sum = J_2_inv_sum + J_2_inv_hat
    

    # Estimate restricted model
    result_cons <- optim(par = param[2], rnd_Probit_LL_cons, y = y, x = x, beta=theta_null,
                         method = c("L-BFGS-B"), lower= c(0),#-Inf,
                         control = list(pgtol=1e-9), hessian=TRUE)
    
    theta_hat_ML_cons = matrix(c(theta_null,result_cons$par), nrow=2)
    test_data[line,"beta_c"]= theta_hat_ML_cons[1]
    test_data[line,"sigma_c"]= theta_hat_ML_cons[2]

  
    ## Score Test
    score = rnd_Probit_LL_g(y,x,theta_hat_ML_cons)
    #J_1 = rnd_Probit_J_1(y,x,theta_hat_ML_cons)
    #print(score)
    #print(J_1)
    #J_1_inv_hat = solve(J_1)/n
    #score = n*t(score)%*% J_1_inv_hat %*%score
    score = n*t(score)%*% J_2_inv_hat %*%score

    test_data[line,"score"]= score
    #J_1_inv_sum = J_1_inv_sum + J_1_inv_hat
    
    
    ## LR Test
    lr = 2*n*(result_cons$value - result$value)
    
    test_data[line,"lr"]= lr
    
  #}
  
  #print(c(J_1_inv_sum,J_2_inv_sum))
  
  return(test_data)
}

##
#n = 300 
#num= 1000
#theta_null <- c(0)
#param = c(0,0)
#cv = qchisq(.95, df=1) 

####      scalar sigmabeta case     |  vector sigmabeta case
###### null   0,0   0,0.05  0,0.15  |  0,0    0,0.05    0,0.15
#beta -2       Y      Y        Y    |   Y       Y          Y         
#beta -1.7                          |                             
#beta -1.5     Y      Y        Y    |   Y       Y          Y      
#beta -1.3                          |                         
#beta -1       Y      Y        Y    |   Y       Y          Y     
#beta -0.7                          |                         
#beta -0.5     Y      Y        Y    |   Y     X_sing     Y completely wrong values            
#beta -0.3                          |                         
#beta 0        Y      Y        Y    | X_sing  X_sing      X_sing           
#beta 0.3                           |                         
#beta 0.5      Y      Y        Y    |   Y       Y          Y     
#beta 0.7                           |                         
#beta 1        Y      Y        Y    |   Y       Y          Y  
#beta 1.3                           |                         
#beta 1.5      Y      Y        Y    |   Y       Y          Y    
#beta 1.7                           |                         
#beta 2        Y      Y        Y    |   Y       Y          Y    
#                                         sigma overestimated

#Dont use params for breakdown
# Also use in-model constraints for testing beta*=0


#a <- probit_beta1_test(n,param,theta_null)
#
#a <- a %>% add_column("reject_wald"=rep(0,num),
#                      "reject_score"=rep(0,num),
#                      "reject_lr"=rep(0,num)) %>% 
#  mutate(reject_wald = if_else(wald > cv,1,0),
#         reject_lr = if_else(lr > cv,1,0))

#mean(a$beta)
#mean(a$sigmasq)

#hist(a$score)

