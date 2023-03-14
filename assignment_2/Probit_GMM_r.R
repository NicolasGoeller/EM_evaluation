Probit_GMM_r <- function(y,x,par,theta1_r) {
  
  par = matrix(c(par,theta1_r),nrow=2, ncol=1)  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  g = (1/n)*(t(y - Phi)%*%x)
  #W_hat = solve((1/n)*t(g)%*%g)
  W_hat = diag(k)
  
  f = (1/n)*(g) %*% W_hat %*% t((1/n)*(g))

  return(f)
}