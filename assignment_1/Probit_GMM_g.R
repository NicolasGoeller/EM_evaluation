Probit_GMM_g <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
	
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  W_hat = diag(k)
  
	f = -2*((1/n)*((t(matrix(rep(phi,k),nrow = n)) * t(x)) %*% x) %*% W_hat %*% ((1/n)*t(x)%*%(y - Phi)))
#                 t(nxk)                        ??t()     nxk    kxk              t(nxk)  nx1
#                         kxk                                       kxk                   kx1
	
	return(f)
}
