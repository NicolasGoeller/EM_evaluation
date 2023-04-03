Probit_LL_r <- function(y,x,par,theta1_r) {
  
  n = length(y) 
  #k = length(par)

  par = matrix(c(par, theta1_r),nrow=2, ncol=1)
  Phi = pnorm(x %*% par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -(1/n)*f
  #f = -f
  return(f)
}
