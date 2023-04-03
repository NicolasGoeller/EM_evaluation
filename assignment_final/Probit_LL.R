Probit_LL <- function(y,x,par) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm(x %*% par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -(1/n)*f

	return(f)
}
