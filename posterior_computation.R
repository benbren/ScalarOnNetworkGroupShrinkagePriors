library(tidyverse)

rinvgamma <- function(n,a,b){
  return(1/rgamma(n,a,b))
}


group_horseshoe_gibs <- function(p, n, burn_ins, posterior_draws = 100, 
                           debug = F,
                           update = "all",
                           show = F,
                           show_params = F,
                           sigma2_beta = 2,
                           prop_zero = 0.8,
                           beta_mean = 5,
                           beta_sd = 0.01,
                           x_cor = 0.5,
                           x_sigma = 1,
                           init_beta = 0
                           ){ 
  
  q = p*(p-1)/2
  dat <- generate_network_data(p = p, n = n,
                               sigma2_beta = sigma2_beta,
                               prop_zero = prop_zero,
                               beta_mean = beta_mean,
                               beta_sd = beta_sd,
                               x_cor = x_cor,
                               x_sigma = x_sigma)
  gam <- data.frame(gam_value = rep(1,p), j = 1:p, a_gam = 1, k = 1:p)
  lambda <- data.frame(lam_value = rep(1, q), j = dat$signal$j, k = dat$signal$k, a_lam = 1)
  sigma2 <- 1
  betas <- data.frame(values  = rep(init_beta, q),  j = dat$signal$j, k = dat$signal$k)
  if(debug){
    betas$values = dat$signal$beta
  }
  tau2 <- 1
  iterations <- burn_ins + posterior_draws
  it <- 0
  beta_mat <- NULL
  sigma_vec <- NULL
  gam_j_vec <- NULL 
  gam_k_vec <- NULL 
  tau_vec <- NULL 
  lam_vec <- 
  X <- NULL 
  for (i in 1:n){
    X <- rbind(X, dat$ntwrks[[i]][lower.tri(dat$ntwrks[[i]])])
  }
 
  
  while(it < iterations){
    
    
    if("atau" %in% update | update == "all"){
      # posterior of a_tau is invgamma(1,1 + 1/tau^2)
      a_tau <- rinvgamma(1,1, 1+1/tau2)
      if(show){
        print(a_tau)
      }
    
    } 
    
    if("agam" %in% update | update == "all"){
      
      for(node in 1:p){
         # posterior of a_gamma_j is invgamma(1, 1 + 1/gamma_j^2)
         gam$a_gam[node]<-rinvgamma(1,1,1 + 1/gam$gam_value[node])
         
      }
    }
    
    if("alam" %in% update | update == "all"){
      
      for(node_pair in 1:q){
        # posterior of a_lam_jk is invgamma(1, 1 + 1/lambda_jk^2)
        lambda$a_lam[node_pair] <- rinvgamma(1,1, 1 + 1/lambda$lam_value[node_pair])
      }
    }
    
    gam_k <- gam %>% select(-j) %>% rename(gam_k_value = gam_value, a_k_gam = a_gam)
    gam_j <- gam %>% select(-k) %>% rename(gam_j_value = gam_value, a_j_gam = a_gam)
    
    
    params <- lambda %>% 
      left_join(gam_k, b = "k") %>% 
      left_join(gam_j, by = "j") %>% 
      left_join(betas, by = c("j","k")) %>% 
      rowwise() %>% 
      mutate(div_sum = values^2/(gam_j_value*gam_k_value*lam_value))
    
    if(show_params){
      print(params)
    }
    

    if("tau" %in% update | update == "all"){
      # posterior of tau^2 is invgamma( (q + 1) / 2, 1 / a_tau + 1/2sigma^2 * sum(beta_jk^2 / lam_jk ^2 gam_j ^2 gam_k ^2))
      tau2 <- rinvgamma(1,(q+1)/2, 1/a_tau + sum(params$div_sum)/(2*sigma2))
      if(show){
        print(tau2)
      }
    }
    
    if("gamma" %in% update | update == "all"){
      for (node in 1:p){
        k_dat <- params %>% filter(k == node | j == node) %>% 
          rowwise() %>% 
          mutate(denom_piece  = ifelse(j == node, (values^2) / (gam_k_value * lam_value) ,values^2 / (gam_j_value * lam_value )))
      
        a_k <- gam %>% filter(k == node) %>% select(a_gam) %>% unlist()
      
        m <- nrow(k_dat)

        s <- sum(k_dat$denom_piece)
      
        # posterior for gamma_k is invgamma(m_k + 1 / 2, 1/a_k + 1/2sigma^2tau^2)
        gam$gam_value[gam$k == node] <- rinvgamma(1,(m+1)/2,1/a_k + s/(2*sigma2*tau2))

      }
    }
    
    if("lambda" %in% update | update == "all"){
    
    for(np in 1:q){
        np_dat <- params[np,]
        j_np <- np_dat$j
        k_np <- np_dat$k
        idx <- lambda$j == j_np & lambda$k == k_np 
      
     
        lambda$lam_value[idx] <- rinvgamma(1,1, 1/np_dat$a_lam + np_dat$values^2/(2*(sigma2*tau2*np_dat$gam_k_value*np_dat$gam_j_value)))
      }
    }
    
    xb <- X%*%betas$values
    
    gam_k <- gam %>% select(-j) %>% rename(gam_k_value = gam_value, a_k_gam = a_gam)
    gam_j <- gam %>% select(-k) %>% rename(gam_j_value = gam_value, a_j_gam = a_gam)
    
    params <- lambda %>% 
      left_join(gam_k, b = "k") %>% 
      left_join(gam_j, by = "j") %>% 
      left_join(betas, by = c("j","k"))
    
    
    LAMBDA_diag <- params %>% rowwise() %>% mutate(diag_el = gam_k_value * gam_j_value*lam_value) %>% select(diag_el) %>% unlist()
    lam <- diag(1,q)
    diag(lam) <- LAMBDA_diag
    inv_lam <- solve(tau2*lam)
    
    
    err = dat$outcomes - xb
    if("sigma" %in% update | update == "all"){
      sigma2 <- rinvgamma(1,(n+q)/2, t(err) %*% err /2 + t(betas$values) %*% (inv_lam) %*% betas$values/2)
      if(show){print(sigma2)}
    }
    
    A <- t(X)%*% X + inv_lam
    
    if(det(A) < 1){
      browser()
    }
  
    Ainv <- solve(A)
    if("betas" %in% update | update == "all"){
      betas$values <- MASS::mvrnorm(1,Ainv%*%t(X)%*%dat$outcomes, sigma2*Ainv)
    }
    print(paste("it:", it))
    
    if(it > burn_ins){
      beta_mat <- cbind(beta_mat, betas$values)
      lam_vec <- cbind(lam_vec, params$lam_value)
      gam_j_vec <- cbind(gam_j_vec, params$gam_j_value)
      gam_k_vec <- cbind(gam_k_vec,params$gam_k_value)
      tau_vec <- c(tau_vec, tau2)
      sigma_vec <- c(sigma_vec, sigma2)
    }

    it <- it + 1
    
    
  }
  
  posterior_pe <- rowMeans(beta_mat)
  
  return(list(real_dat = dat, posterior_draws_beta = beta_mat, 
              posterior_pe = posterior_pe, 
         posterior_draws_lam = lam_vec, 
         posterior_draws_gam_j = gam_j_vec,
         posterior_draws_gam_k = gam_k_vec, 
         posterior_draws_tau = tau_vec, 
         posterior_draws_sigma = sigma_vec))
}

