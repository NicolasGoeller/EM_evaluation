probit_beta_test <- function(n,num,param,theta_null){
  
  set.seed(1) #Seed for random number generator
  
  ## Load functions
  source("assignment_final/rnd_Probit_LL.R")
  source("assignment_final/rnd_Probit_LL_cons.R")
  source("assignment_final/rnd_Probit_LL_g.R")
  source("assignment_final/rnd_Probit_J_1.R")
 
  test_data <- data.frame("beta_star"=rep(param[1],num),"sigma_star"=rep(param[2],num),
                          "wald"=rep(NA,num),"score"=rep(NA,num),"lr"=rep(NA,num),
                          "beta"=rep(NA,num),"sigma"=rep(NA,num),
                          "beta_c"=rep(NA,num),"sigma_c"=rep(NA,num))

  k = length(param)
  i = 1
  
  #J_1_inv_sum = matrix(c(0,0,0,0),ncol=k)
  J_2_inv_sum = matrix(c(0,0,0,0),ncol=k)
  
  while (i <= num) {

    ## Random Coefficient Probit DGP
    x = matrix(rnorm((k-1)*n,0,1),ncol=k-1) 	# regressors

    u = rnorm(n,0,1)                            # error term 1
    e = rnorm(n,0,1)                            # error term 2
    beta_i = param[1]+param[2]*e                # varying coefficients
    y_star = x*beta_i + u                       # latent "utility"

    y = ceiling(y_star/(max(abs(y_star))+0.1))

    ## Estimate unrestricted model
    result <- tryCatch({
      optim(par = param, rnd_Probit_LL, y = y, x = x,
            method = c("L-BFGS-B"), lower= c(-Inf,0),
            control = list(pgtol=1e-9), hessian=TRUE)},
      error= function(e){
        return(NULL)
      })
    if (is.null(result) == TRUE) {
      print("Repeat")
      next
    }

    theta_hat_ML = result$par
    test_data[i,"beta"]= theta_hat_ML[1]
    test_data[i,"sigma"]= theta_hat_ML[2]
    

    ## Wald Test
    J_2_inv_hat = tryCatch({
      solve(result$hessian)/n},
      error = function(e){
        return(NULL)
    })
    if (is.null(J_2_inv_hat) == TRUE) {
      print("Repeat")
      next
    }
    
    wald = (theta_hat_ML[1] - theta_null)^2/J_2_inv_hat[1,1]
    
    test_data[i,"wald"]= wald
    J_2_inv_sum = J_2_inv_sum + J_2_inv_hat
    
    
    # Estimate restricted model
    result_cons <- tryCatch({
      optim(par = param[2], rnd_Probit_LL_cons, y = y, x = x, beta=theta_null,
            method = c("L-BFGS-B"), lower= c(0),#-Inf,
            control = list(pgtol=1e-9), hessian=TRUE)},
      error = function(e){
        return(NULL)
        })
    if (is.null(result_cons) == TRUE) {
      print("Repeat")
      next
    } 
    
    theta_hat_ML_cons = matrix(c(theta_null,result_cons$par), nrow=2)
    test_data[i,"beta_c"]= theta_hat_ML_cons[1]
    test_data[i,"sigma_c"]= theta_hat_ML_cons[2]

  
    ## Score Test
    score = rnd_Probit_LL_g(y,x,theta_hat_ML_cons)
    #J_1 = rnd_Probit_J_1(y,x,theta_hat_ML_cons)
    #J_1_inv_hat = solve(J_1)/n
    #score = n*t(score)%*% J_1_inv_hat %*%score
    score = n*t(score)%*% J_2_inv_hat %*%score

    test_data[i,"score"]= score
    #J_1_inv_sum = J_1_inv_sum + J_1_inv_hat
    
    
    ## LR Test
    lr = 2*n*(result_cons$value - result$value)
    
    test_data[i,"lr"]= lr
    
    i = i+1
    
  }
  
  #print(c(J_1_inv_sum,J_2_inv_sum))
  
  return(test_data)
}

##
#n = 300 
#num= 1000
#theta_null <- c(0)
#param = c(0,0.05)
#cv = qchisq(.95, df=1) 

#tic()
#a <- probit_beta1_test(n,num,param,theta_null)
#toc()
#
#a <- a %>% add_column("reject_wald"=rep(0,num),
#                      "reject_score"=rep(0,num),
#                      "reject_lr"=rep(0,num)) %>% 
#  mutate(reject_wald = if_else(wald > cv,1,0),
#         reject_lr = if_else(lr > cv,1,0))

#mean(a$beta)
#mean(a$sigmasq)

#hist(a$score)

