rnd_Probit_LL_g <- function (y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  sigmabeta = sqrt(1 + (x^2)*par[2]^2)
  sigmabeta2 = sqrt(1 + (x^2)*par[2]^2)^1.5
  
  Phi = pnorm((x * par[1])/sigmabeta) # generate cdf
  phi = dnorm((x * par[1])/sigmabeta) # generate pdf
  
  d_beta = rep(NA,n)
  d_sigma = rep(NA,n)
  
  for (i in 1:n){
    d_beta[i] = (y[i]*x[i]*phi[i])/(sigmabeta[i]*Phi[i]) - ((1-y[i])*x[i]*phi[i])/((1-Phi[i])*sigmabeta[i])
    d_sigma[i] = ((1-y[i])*par[1]*x[i]^3*phi[i])/(2*(1-Phi[i])*sigmabeta2[i]) - (y[i]*par[1]*x[i]^3*phi[i])/(2*Phi[i]*sigmabeta2[i])
  }

  g = matrix(c(sum(d_beta),sum(d_sigma)),nrow=k)
  g = -(1/n)*g
  
	return(g)
}