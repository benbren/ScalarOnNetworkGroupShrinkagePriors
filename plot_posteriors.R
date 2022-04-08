summarize_posteriors <- function(sim){
  signals <- sim$real_dat$signal 
  p <- max(signals$j)
  draws <- dim(sim$posterior_draws_beta)[2]
  beta_g <- cbind(sim$posterior_draws_beta, signals) %>% 
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,p,p,0)) %>% 
            as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j ~ k) + 
    theme_bw() + xlab(" Posterior Draws (Beta)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") 
  
  lam_g <- cbind(sim$posterior_draws_lam, signals) %>%
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,p,p,0)) %>% 
    as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j ~ k) + 
    theme_bw() + xlab(" Posterior Draws (Edge Shrinkage)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") + xlim(0,25)
  
  gam_g <- cbind(sim$posterior_draws_gam_j, signals) %>%
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,p,p,0)) %>% 
    as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    mutate(edge = paste0(j,", ", k)) %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j~k) + 
    theme_bw() + xlab(" Posterior Draws (Node Shrinkage)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") 
  
  post_lam <- cbind(sim$posterior_draws_lam,signals) %>% as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    group_by(j,k) %>% summarize(pe_lam_jk = mean(posterior_draws))
  
  post_gam_j <- cbind(sim$posterior_draws_gam_j,signals) %>% as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    group_by(j,k) %>% summarize(pe_gam_j = mean(posterior_draws)) 
  
  post_gam_k <- cbind(sim$posterior_draws_gam_k,signals) %>% as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    group_by(j,k) %>% summarize(pe_gam_k = mean(posterior_draws)) 
  
  post_gam <- full_join(post_gam_j, post_gam_k, by = c("j","k")) 
  
  tau_plot <- hist(sim$posterior_draws_tau)
  sigma_plot <- hist(sim$posterior_draws_sigma)
  
  return(list("beta_plot" = beta_g, "lam_plot" = lam_g, "gam_plot" = gam_g,
               "tau_plot" = tau_plot, "sigma_plot" = sigma_plot,
              "lam_df"= post_lam, 
              "gam_df" = post_gam,
              "post_sigma" = mean(sim$posterior_draws_sigma),
              "post_tau" = mean(sim$posterior_draws_tau)))
  
} 



# also want to compare with emily's method, and linear regression, ridge? and LASSO. 
#look for software of other new methods
# MSE for betas, selection accuracy, squared prediction accuracy  

# wnt to consdier vriation in simulation ie sample size (small = 100), variations in signal to not ratio
# var(xb) / var(Y) - large means a stronng signal 
# and number of nodes (15 = small, large = small = make very sparse)  



# we want to extend to logistic regression (can do from paper on sampling)
