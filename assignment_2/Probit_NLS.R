Probit_NLS <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
	
  Phi = pnorm(x %*% par)
  
  f = sum((y - Phi)^2)
  f = f/n

	return(f)
}
