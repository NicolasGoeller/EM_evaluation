Probit_NLS_g <- function (y,x,par) {

	n = length(y) 
	k = length(par)
	
	Phi = pnorm(x %*% par) # generate cdf
	phi = dnorm(x %*% par) # generate pdf

	f = t((y-Phi)*phi) %*% x
	f = -(2/n)*f

	return(f)
}