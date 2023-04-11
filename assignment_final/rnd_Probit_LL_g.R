rnd_Probit_LL_g <- function (y,x,par) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm((x %*% par[1:2])/sqrt(1 + (x[,2]**2)*(par[3]**2))) # generate cdf
  phi = dnorm((x %*% par[1:2])/sqrt(1 + (x[,2]**2)*(par[3]**2))) # generate pdf
  
  
  f = t(y*phi/Phi - (1-y)*phi/(1-Phi)) %*% x
  f = -(1/n)*t(f)
  
	return(f)
}