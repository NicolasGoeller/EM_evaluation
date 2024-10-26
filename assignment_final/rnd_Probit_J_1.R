rnd_Probit_J_1 <- function(y,x,par) {

  n = length(y)
  k = length(par)
  
  sigmabeta = sqrt(1 + (x^2)*par[2]^2)
  sigmabeta2 = sqrt(1 + (x^2)*par[2]^2)^1.5

  Phi = pnorm((x * par[1])/sigmabeta[1]) # generate cdf
  phi = dnorm((x * par[1])/sigmabeta[1]) # generate pdf
  
  d_beta = rep(NA,n)
  d_sigma = rep(NA,n)

  for (i in 1:n){
    d_beta[i] = (y[i]*x[i]*phi[i])/(sigmabeta[1]*Phi[i]) - ((1-y[i])*x[i]*phi[i])/((1-Phi[i])*sigmabeta[1])
    d_sigma[i] = ((1-y[i])*par[1]*x[i]^3*phi[i])/(2*(1-Phi[i])*sigmabeta2[1]) - (y[i]*par[1]*x[i]^3*phi[i])/(2*Phi[i]*sigmabeta2[1])
    #print(c(d_beta[i], d_sigma[i]))
  }

  score = matrix(rbind(d_beta,d_sigma), nrow=k)
  
  f = score%*%t(score)

	return(f)
}
