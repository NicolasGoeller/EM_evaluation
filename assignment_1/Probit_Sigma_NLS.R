Probit_Sigma_NLS <- function(y,x,par) {

  ## According to wooldridge, \Sigma = Var(s(w,\theta)) variance of score function
  ## Score function is the expectation of the gradient of objective function (by theta)
  
  ## H is the expectation of the Hessian matrix for true theta
  
    
  n = length(y)
  k = length(par)
  
  Phi = pnorm(x %*% par) # generate cdf
  phi = dnorm(x %*% par) # generate pdf
  
  #s = t((y-Phi)*phi) %*% x
  s = t(x) %*% ((y-Phi)*phi) 
  s = -(2/n)*s
  
  
  
  #expectation of square of score function
  #f = (1/n)*(t(s) %*% s) #this gives 2x2 matrix
  f = (1/n)*(s %*% t(s)) #this gives 2x2 matrix

	return(f)
}
