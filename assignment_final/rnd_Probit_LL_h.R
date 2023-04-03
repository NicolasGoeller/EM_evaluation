Probit_LL_h <- function (y,x,par) {
  
	Phi = pnorm(x %*% par)
	phi = dnorm(x %*% par)

	n = length(y) 
	k = length(par)

	# Computing the Hessian
	h = t(x)%*% (matrix(rep((1-y)*(phi/(1-Phi)*(x%*%par)-(phi^2)/(1-Phi)^2),k),nrow=n)*x) + 
		t(x)%*% (matrix(rep((y)*(phi/Phi*(-x%*%par)-(phi^2)/Phi^2),k),nrow=n)*x)
	h = -h

	return(h)
}