mcprobit_beta1_test <- function(n,num,param){
  
  set.seed(1) #Seed for random number generator
  
  ## Load functions
  source("assignment_final/rnd_Probit_LL.R")
  source("assignment_final/rnd_Probit_LL_g.R")
  source("assignment_final/rnd_Probit_J_1.R")
 
  theta0 = list(mean= c(1,param[1]),
                var= c(0,param[2])) #Random coefficient model: beta0, beta1, sigma
  
  test_data <- data.frame("beta1"=rep(param[1],num),"sigma"=rep(param[2],num),
                          "wald"=rep(NA,num),"score"=rep(NA,num),"lr"=rep(NA,num))
  
  for (i in 1:num){
    ## Random Coefficient Probit DGP
    x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
    u = rnorm(n,0,1 + (x**2)%*%(theta0[[2]]**2))				# error term
    y_star = x %*% theta0[[1]] + u						# latent "utility"
    y = ceiling(y_star/(max(abs(y_star))+0.1))
    
    ## Wald Test
    wald = (theta_hat_ML[2] - theta1_null)^2/J_2_inv_hat[2,2]
    
    test_data[i,"wald"]= wald
    
    ## Score Test
    score = Probit_LL_g(y,x,theta_hat_ML_cons)
    score = n*score%*%solve(Probit_J_1(y,x,theta_hat_ML_cons)/n)%*%t(score)
    
    test_data[i,"score"]= score
    
    ## LR Test
    lr = 2*n*(result_cons$value - result$value)
    
    test_data[i,"lr"]= lr

    
  }
  return(test_data)
}