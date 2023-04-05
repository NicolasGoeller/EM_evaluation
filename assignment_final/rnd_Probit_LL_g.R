Probit_LL_g <- function (y,x,par1,par2) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm((x %*% par1)/sqrt(1 + (x**2)%*%(par2**2))) # generate cdf
  phi = dnorm((x %*% par1)/sqrt(1 + (x**2)%*%(par2**2))) # generate pdf
  
  
  f = t(y*phi/Phi - (1-y)*phi/(1-Phi)) %*% x
  f = -(1/n)*t(f)
  
	return(f)
}