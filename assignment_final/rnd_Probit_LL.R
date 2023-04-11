rnd_Probit_LL <- function(y,x,par) {
  
  n = length(y) 
  #k = length(par)
  
  Phi = pnorm((x %*% par[1:2])/sqrt(1 + (x[,2]**2)*(par[3]**2)))
  #print(sum(is.na(Phi)))
  #print(dim(Phi))
  print(c(min(Phi),mean(Phi),max(Phi)))
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  #print(sum((1-y)*log(1-Phi)))
  #print(f)
  f = -(1/n)*f
  
  if (is.na(f)){
    print(f)
    }
  
	return(f)
}

