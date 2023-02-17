Probit_Sigma_NLS <- function(y,x,par) {

  ## According to wooldridge, \Sigma = Var(s(w,\theta)) variance of score function
  ## Score function is the expectation of the gradient of objective function (by theta)
  
    
  n = length(y)
  k = length(par)
  
  ...
  
  f = ...

	return(f)
}
