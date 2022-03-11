generate_network_data <- function(p, 
                                  n,
                                  sigma2_beta = 2,
                                  prop_zero = 0.8,
                                  beta_mean = 5,
                                  beta_sd = 0.01,
                                  x_cor = 0.5,
                                  x_sigma = 1){
  
  
  ## Create a random Adjacency Matrix that represents the network 
  adjm <- matrix(NA,p,p)
  q = p*(p-1)/2
  true_signals <- matrix(sample(0:1,q , replace = T, prob = c(prop_zero,1-prop_zero)))
  adjm[lower.tri(adjm)] <- true_signals
  adjm[upper.tri(adjm)] <- t(true_signals)
  diag(adjm) <- 0
  
  ntwrk <- igraph::graph_from_adjacency_matrix(adjm,
                                               mode= "undirected")
  

  
  i <- 2
  j <- NULL # rows 
  k <- NULL # columns 
  
  while (i < p+1){
    addonj <- i:p
    addonk <- rep(i-1,p+1-i) 
    j <- c(j,addonj)
    k <- c(k,addonk)
    i <- i + 1
  }
  
  beta_with_indices = data.frame(signal = true_signals, j = j, k = k)
  beta_with_indices$beta <- beta_with_indices$signal * rnorm(dim(beta_with_indices)[1],beta_mean, beta_sd)
  # 
  

  cov_mat <- diag(q)*x_sigma
  cov_mat[upper.tri(cov_mat)] <- x_cor
  cov_mat[lower.tri(cov_mat)] <- x_cor
  observed_networks <- list()
  observed_networkslt <- list()
  
  y <- NULL
  
  for (i in 1:n){
    x_i <- MASS::mvrnorm(n = 1, mu = rep(0,q), Sigma = cov_mat)
    ob_ntwrk <- matrix(NA,p,p)
    ob_ntwrk[upper.tri(ob_ntwrk)] <- x_i
    ob_ntwrk[lower.tri(ob_ntwrk)] <- t(x_i)
    diag(ob_ntwrk) <- 0
    
    y_i <- t(beta_with_indices$beta)%*%ob_ntwrk[lower.tri(ob_ntwrk)] + rnorm(1,mean = 0, sd = sqrt(sigma2_beta))
    # 
    y <- c(y, y_i)
    
    observed_networkslt[[i]] <- ob_ntwrk[lower.tri(ob_ntwrk)] # make sure this matches indices!! 
    observed_networks[[i]] <- ob_ntwrk # make sure this matches indices!! 
  }
  
  rtrn <- list(outcomes = y, signal = beta_with_indices, ntwrks = observed_networks, ntworkslt = observed_networkslt)
  return(rtrn)
  
  # memory efficient is to save onky upper triangular 
  
  # to consider: implement posterior, only need to work on upper triangular 
  # carefully consider the indices with the priors 
}

  

