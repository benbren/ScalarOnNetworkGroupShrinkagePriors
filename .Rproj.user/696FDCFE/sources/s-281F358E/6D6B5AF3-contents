generate_network_data <- function(p, 
                                  n,
                                  sigma2_beta = 2,
                                  prop_zero = 0.8,
                                  beta_mean = 5,
                                  beta_sd = 0.01,
                                  x_cor = 0.5,
                                  x_sigma = 1,
                                  beta_hubs = F,
                                  hub_nodes = NULL,
                                  hub_degrees = NULL){
  
  #profvis({
  ## Create a random Adjacency Matrix that represents the network 
  q = p*(p-1)/2
  if(!beta_hubs){
    adjm <- matrix(NA,p,p)
    true_signals <- matrix(sample(0:1,q , replace = T, prob = c(prop_zero,1-prop_zero)))
    adjm[lower.tri(adjm)] <- true_signals
    adjm[upper.tri(adjm)] <- t(true_signals)
    diag(adjm) <- 0
    
    ntwrk <- igraph::graph_from_adjacency_matrix(adjm,
                                                 mode= "undirected")
  } else {
    if(is.null(hub_nodes | is.null(hub_degrees))){
      error("Please specify nodes AND degress for the hubs")
    }
    hub_connects <- list()
    if(p - length(hub_nodes) - sum(hub_degrees) < 2){
      error("Pick a larger p - your hub and connections are greater than the number of nodes ")
    }
    num_nodes <- p
    adjm <- matrix(0,nrow=num_nodes, ncol=num_nodes)
    for(i in 1:length(hub_nodes)){
      if(i == 1)
        hub_connects[[i]] <- length(hub_nodes)+1:hub_degrees[1]
      else
        hub_connects[[i]] <- max(hub_connects[[i-1]]) + 1:hub_degrees[i]
      adjm[hub_nodes[i],hub_connects[[i]]] <- 1
    }
    adjm <- t(adjm) + adjm
    
    true_signals <-  adjm[lower.tri(adjm)]
    diag(adjm) <- 0 
    
    ntwrk <- igraph::graph_from_adjacency_matrix(adjm,mode="undirected")
  }
  
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


  #cov_mat <- diag(q)*x_sigma
  cov_mat <- toeplitz(c(x_sigma, rep(x_cor,q-1)))
  #cov_mat[upper.tri(cov_mat)] <- x_cor
  #cov_mat[lower.tri(cov_mat)] <- x_cor
  observed_networks <- list()
  observed_networkslt <- list()
  
  y <- NULL
  x <- mvnfast::rmvn(n,rep(0,q), cov_mat, ncores = 2)
  for (i in 1:n){
    x_i <- x[i,]
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
  
  rtrn <- list(outcomes = y, signal = beta_with_indices, ntworkslt = observed_networkslt,
               p = p, n = n, full_ntwork = ntwrk)
  #})
  return(rtrn)
  
}

  

