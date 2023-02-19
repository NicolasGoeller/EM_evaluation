Probit_Var_GMM <- function(y,x,par) {
  
  ## G_hat is the gradient of g(z,theta) at theta_hat
  ## Omega_hat is the variance of the g(z,theta) function - sample mean of g()g()'
  
  n = length(y)
  k = length(par)

  Phi = pnorm(x %*% par) # generate cdf
  phi = dnorm(x %*% par) # generate pdf
  W_hat = diag(k)
  
  s = t(y-Phi) %*% x
  Omega_hat = (1/n)*(t(s) %*% s)
    
  G_hat = (1/n) * (t(matrix(rep(-phi,k),nrow=n)*x) %*% x)
  
  f = solve(G_hat) %*% Omega_hat %*% solve(t(G_hat))
  
	return(f)
}

