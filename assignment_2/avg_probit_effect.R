avg_probit_effect <- function(theta_vec, theta_pos, x){
  
  n = dim(x)[1]
  phi = dnorm(x %*% theta_vec)
  
  avg_effect = sum(phi*theta_vec[theta_pos+1])
  avg_effect = (1/n)*avg_effect
  
  return(avg_effect)
}