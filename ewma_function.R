ewma <- function(rtn = rtn, lambda = .6, increment = .000000001, list_out = TRUE) {

  vcov_list <- vector(mode = "list", length = (dim(rtn)[1]))
  if (is.null(row.names(rtn)) == FALSE) names(vcov_list) <- row.names(rtn)

  vcov_init <- cov(rtn)
  decomp <- eigen(vcov_init)
  offset <- 0
  while (base::min(decomp$values) < 0) {
    offset <- offset + increment
    V <- diag(decomp$values) + offset * diag(dim(vcov_init)[1])
    vcov_init <- decomp$vectors %*% V %*% t(decomp$vectors)
    decomp <- eigen(vcov_init)
  }
  print(paste0("offset used on diagonal elements of initial cov matrix to ensure PD: ", offset))
  
  for (t in 1:dim(rtn)[1]) {
    rtt <- matrix(rep(rtn[t, ], dim(rtn)[2]), nrow = dim(rtn)[2], byrow = TRUE)
    if (t == 1) {
      vcov_list[[t]] <- (1 - lambda) * rtt * rtn[t, ] + lambda * vcov_init 
      vcov_mat <- c(vcov_list[[t]])
    } else {
      vcov_list[[t]] <- (1 - lambda) * rtt * rtn[t, ] + lambda * vcov_list[[t - 1]]
      vcov_mat <- rbind(vcov_mat, c(vcov_list[[t]]))
    } 
    if (is.null(colnames(rtn)) == FALSE) {
      colnames(vcov_list[[t]]) <- colnames(rtn)
      row.names(vcov_list[[t]]) <- colnames(rtn)
    }   
  }

  if (list_out == FALSE) return(vcov_mat) else return(vcov_list)
  
}

  