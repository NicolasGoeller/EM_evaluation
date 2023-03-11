Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
	
  Phi = pnorm(x %*% par)
  W_hat = diag(k)
  
	f = (1/n)*(t(y - Phi)%*%x) %*% W_hat %*% t((1/n)*(t(y - Phi)%*%x))

	return(f)
}
