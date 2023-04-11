mcprobit_beta1_test <- function(n,num,param,theta_null){
  
  set.seed(1) #Seed for random number generator
  
  ## Load functions
  source("assignment_final/rnd_Probit_LL.R")
  source("assignment_final/rnd_Probit_LL_g.R")
  source("assignment_final/rnd_Probit_J_1.R")
 
  test_data <- data.frame("beta1"=rep(param[2],num),"sigma"=rep(param[3],num),
                          "wald"=rep(NA,num),"score"=rep(NA,num),"lr"=rep(NA,num))
  
  k = length(param)-1

  
  for (i in 1:num){
    print(i)
    ## Random Coefficient Probit DGP
    x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
    u = rnorm(n,0,1)
    
    #e = rnorm(n,0,1)
    #beta_i = param[2] + x[,2]*e
    #y_star = param[1] + beta_1*x[,2] + u
    
    u = rnorm(n,0,1 + (x[,2]**2)*(param[3]**2))				# error term
    y_star = x %*% param[1:2] + u						# latent "utility"
    y = ceiling(y_star/(max(abs(y_star))+0.1))
    
    ## Estimate unrestricted model
    result <- optim(par = theta_null, rnd_Probit_LL, y = y, x = x, 
                    method = c("L-BFGS-B"), lower= c(-Inf,-Inf,0),
                    control = list(pgtol=1e-9), hessian=TRUE)

    theta_hat_ML = result$par
    
    ## Wald Test
    #wald = (theta_hat_ML[2] - theta1_null)^2/J_2_inv_hat[2,2]
    wald = t(theta_hat_ML - theta_null) %*% result$hessian %*% (theta_hat_ML - theta_null)
    
    test_data[i,"wald"]= wald
    
    
    # Estimate restricted model
    result_cons <- optim(par = theta_null, rnd_Probit_LL, y = y, x = x,
                         method = c("L-BFGS-B"), lower= c(-Inf,0,0),
                         control = list(pgtol=1e-9), hessian=TRUE)
    theta_hat_ML_cons = result_cons$par
    
    
    ## Score Test
    #score = rnd_Probit_LL_g(y,x,theta_hat_ML_cons)
    #print("test4")
    #score = n*score%*%solve(rnd_Probit_J_1(y,x,theta_hat_ML_cons)/n)%*%t(score)
    
    test_data[i,"score"]= NA#score
    
    ## LR Test
    lr = 2*n*(result_cons$value - result$value)
    
    test_data[i,"lr"]= lr

    
  }
  return(test_data)
}

##
n = 300 
num= 100
theta_null <- c(1,0,0.05)
param = c(1,1.3,0) 

mcprobit_beta1_test(n,num,param,theta_null)
