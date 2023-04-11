rnd_Probit_J_1 <- function(y,x,par) {
  
  n = length(y)
  m = length(par)
  
  Phi = pnorm((x %*% par[1:2])/sqrt(1 + (x[,2]**2)*(par[3]**2)))
  phi = dnorm((x %*% par[1:2])/sqrt(1 + (x[,2]**2)*(par[3]**2)))
  
  print(dim(phi))
  print(dim(Phi))
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),m),nrow=n)*x
  
  f = (1/n)*t(g)%*%g

	return(f)
}
