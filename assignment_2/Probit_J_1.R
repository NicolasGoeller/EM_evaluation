Probit_J_1 <- function(y,x,par) {
  
  n = length(y)
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = (1/n)*t(g)%*%g

	return(f)
}
