rnd_Probit_LL <- function(y,x,par) {
  
  n = length(y) 
  
  Phi = pnorm((x * par[1])/sqrt(1 + (x^2)*par[2]^2))
  
  ## Normalise edges of logged Phi vectors
  log1 = log(Phi)
  log1[log1 == -Inf] = -400
  log1[log1 == Inf] = 400
  
  log2 = log(1-Phi)
  log2[log2 == -Inf] = -400
  log2[log2 == Inf] = 400
  
  # Computing the log-likelihood
  f = sum(y*log1) + sum((1-y)*log2)

  f = -(1/n)*f
  #print(paste(par[1], par[2], f))
  
	return(f)
}



