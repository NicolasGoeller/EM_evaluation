Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
	
  Phi = pnorm(x %*% par)
  W_hat = diag(k)
  
  print(dim((1/n)*(t(y - Phi)%*%x)))
  print(dim(W_hat))
  print(dim(t((1/n)*(t(y - Phi)%*%x))))
  
	f = (1/n)*(t(y - Phi)%*%x) %*% W_hat %*% t((1/n)*(t(y - Phi)%*%x))
	f = -f

	return(f)
}
