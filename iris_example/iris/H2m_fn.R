
# This is the density function for the integral
# rearranges matrix inverses and uses Student-t CDF
H2m_int_cdf <- function(S, beta, ds) {
  m <- nrow(ds)
  nu <- 2*beta-1
  s <- qt(S, nu)/sqrt(nu)
  h_val <- 1
  for (ii in 1:m) {
    quad_part <- ds[ii,1]^2 + ds[ii,2]^2
    local_part <- (1 - (quad_part - 2*1i*ds[ii,1]*ds[ii,2]*s)/(1+s^2))/(1-quad_part)
    h_val <- h_val*local_part^(-beta+(m-1)/2)
  }
  Re(h_val)
}

# This is the H_{2,m} integral
# Main arguments are the parameter beta and ds, the matrix of rows from the data
# by default returns the log value, and passes tolerance tol to the numerical integration
H2m_fn <- function(beta, ds, log = TRUE, tol = 1e-12) {
  m <- nrow(ds) # how many determinants
  overall_factor <- (m*(m-1)/4+1/2)*log(pi) + m*lgamma(beta - (m-1)/2) - 2*lgamma(beta) - (beta-(m-1)/2)*sum(log(1-rowSums(ds^2)))
  if (m>2) {
    overall_factor <- overall_factor - sum(lgamma(beta - (3:m-1)/2))
  }
  H_result <- log(integrate(H2m_int_cdf, lower = 0, upper = 1, beta = beta, ds = ds, rel.tol = tol)$value) 
  if (log) {
    return(H_result + overall_factor)
  } else {
    return(exp(H_result + overall_factor))
  }
}

