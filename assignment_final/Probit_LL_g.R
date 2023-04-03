Probit_LL_g <- function (y,x,par) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm(x %*% par) # generate cdf
  phi = dnorm(x %*% par) # generate pdf
  
  
  f = t(y*phi/Phi - (1-y)*phi/(1-Phi)) %*% x
  f = -(1/n)*t(f)
  
	return(f)
}