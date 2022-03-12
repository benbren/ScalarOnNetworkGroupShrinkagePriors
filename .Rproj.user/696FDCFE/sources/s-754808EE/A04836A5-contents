plot_posteriors <- function(simulation){
  signals <- sim$real_dat$signal 
  draws <- dim(sim$posterior_draws_beta)[2]
  beta_g <- cbind(sim$posterior_draws_beta, signals) %>% 
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,7,7,0)) %>% 
            as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j ~ k) + 
    theme_bw() + xlab(" Posterior Draws (Beta)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") 
  
  lam_g <- cbind(sim$posterior_draws_lam, signals) %>%
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,7,7,0)) %>% 
    as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j ~ k) + 
    theme_bw() + xlab(" Posterior Draws (Edge Shrinkage)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") + xlim(0,25)
  
  gam_g <- cbind(sim$posterior_draws_gam_j, signals) %>%
    rbind(c(rep(NA,draws),0,1,1,0)) %>%
    rbind(c(rep(NA,draws),0,7,7,0)) %>% 
    as.tbl() %>% 
    pivot_longer(-c(signal,j,k,beta), values_to = "posterior_draws") %>% 
    mutate(edge = paste0(j,", ", k)) %>% 
    ggplot(aes(x = posterior_draws, color = as.factor(signal))) + geom_density() + facet_grid(j~k) + 
    theme_bw() + xlab(" Posterior Draws (Node Shrinkage)") + ylab("Density") + 
    theme(legend.box = "horizontal") + labs(color = "True Signal") 
  
  tau_plot <- hist(sim$posterior_draws_tau)
  sigma_plot <- hist(sim$posterior_draws_sigma)
  
  return(list("beta_plot" = beta_g, "lam_plot" = lam_g, "gam_plot" = gam_g,
               "tau_plot" = tau_plot, "sigma_plot" = sigma_plot))
  
} 
