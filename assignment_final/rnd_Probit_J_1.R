Probit_J_1 <- function(y,x,par1,par2) {
  
  n = length(y)
  k = length(par)
  
  Phi = pnorm((x %*% par1)/sqrt(1 + (x**2)%*%(par2**2)))
  phi = dnorm((x %*% par1)/sqrt(1 + (x**2)%*%(par2**2)))
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = (1/n)*t(g)%*%g#

	return(f)
}
