Probit_Var_GMM <- function(y,x,par) {
  
  n = length(y)
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  m = matrix(rep((y-Phi),k),nrow=n)*x
  A = t(x) %*% (matrix(rep(phi,k),nrow=n)*x)
  A_inv = solve(A)
  
  var_hat = A_inv %*% (t(m) %*% m) %*% A_inv
  
  return(var_hat)
}