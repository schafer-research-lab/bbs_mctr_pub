ewma_weights <- function(lambda = lambda_val, mval = 50, weight_cutoff = .993) {
  ew <- ewc <- NA
  for (t in 1:mval) {
    ew[t] <- (1 - lambda) * lambda^(mval - t) 
  }  
  ewc <- cumsum(rev(ew))  
  num_years <- base::min(which(ewc > weight_cutoff))
  return(num_years)
}
