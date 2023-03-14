GMM_QLR <- function(y,x,par) {
  
  n = length(y)
  
  Phi = pnorm(x %*% par)
  
  g = (1/n)*(t((y-Phi))%*%x)
  
  Omega_hat = GMM_Omega(y,x,par)

  QLR_hat = n*g%*%solve(Omega_hat)%*%t(g)
  
  return(QLR_hat)
}
