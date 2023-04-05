Probit_LL <- function(y,x,par1,par2) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm((x %*% par1)/sqrt(1 + (x**2)%*%(par2**2)))
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -(1/n)*f

	return(f)
}
