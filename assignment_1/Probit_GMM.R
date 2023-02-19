Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
	
  Phi = pnorm(x %*% par)
  W_hat = diag(k)
  
	f = (1/n)*(t(y - Phi)%*%x) %*% W_hat %*% t((1/n)*(t(y - Phi)%*%x))
	#f = (1/n^2)*((t(y - Phi)%*%x) %*% W_hat %*% (t(x) %*% (y - Phi))) #Andrea's version
	#f = -f

	return(f)
}
