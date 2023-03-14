GMM_Omega <- function(y,x,par){
  
  n = length(y)
  
  Phi = pnorm(x %*% par)
  
  Omega_diag1 = (1/n)*sum((y-Phi)^2)
  Omega2_offdiag = (1/n)*sum((y-Phi)^2 * x[,2])
  Omega_diag2 = (1/n)*sum((y-Phi)^2 * x[,2]^2)
  Omega_hat = matrix(data=c(c(Omega_diag1,Omega2_offdiag),
                            c(Omega2_offdiag,Omega_diag2)), ncol=2)
  
  return(Omega_hat)
}