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

plot_gammas <- function(ghs_dat){
  
  gam_dat <- data.frame(gam_j_mean = rowMeans(ghs_dat$posterior_draws_gam_j),
                        gam_k_mean = rowMeans(ghs_dat$posterior_draws_gam_k),
                        j = ghs_dat$real_dat$signal$j, 
                        k = ghs_dat$real_dat$signal$k,
                        sig = ghs_dat$real_dat$signal$signal)
  
  gam_dat1 <- gam_dat %>% pivot_longer(-c(gam_k_mean, gam_j_mean, sig)) 
  connections <- gam_dat1 %>% group_by(value) %>% summarise(connects = sum(sig), sig = ifelse(any(sig == 1),TRUE,FALSE))
  plt <- gam_dat1 %>% select(-sig) %>% distinct(value, .keep_all = T) %>% 
    arrange(value) %>% mutate(gam_j_mean = ifelse(name == "k", gam_k_mean, gam_j_mean)) %>%
    select(-name, - gam_k_mean) %>% left_join(connections) %>% 
    rename(node = value, gam_mean = gam_j_mean) %>% 
    ggplot(aes(x = node, y = gam_mean, color = as.factor(sig), shape = as.factor(connects))) + 
    geom_point(size = 4) + theme_bw() + 
    labs(shape = "Node Degrees", y = "Mean Posterior Node-Specific Shrinkage Parameter", 
         x = "Node", color = "True Signal") + 
    theme(text=element_text(size=7.5))
  
  return(plt)
  
}

plot_lams <- function(ghs_dat){
  
  lam_dat <- data.frame(lam_mean = rowMeans(ghs_dat$posterior_draws_lam),
                        j = ghs_dat$real_dat$signal$j, 
                        k = ghs_dat$real_dat$signal$k,
                        sig = ghs_dat$real_dat$signal$signal
  )
  
  plt <- lam_dat %>% ggplot(aes(x = as.factor(k), y = as.factor(j), color = as.factor(sig), fill = as.factor(sig), alpha = lam_mean)) + 
    geom_tile() + theme_bw() + labs(alpha = "Mean Posterior Edge-Specific Shrinkage",
                                    fill = "True Signal", 
                                    x = "Row Node",
                                    y = "Column Node") + guides(color = FALSE) + 
    theme(text=element_text(size=7.5))
  
  return(plt)
}

plot_tau <- function(ghs_dat){
  
  tau_dat <- data.frame(tau = ghs_dat$posterior_draws_tau, 
                        id = 1:length(ghs_dat$posterior_draws_tau)) 
  plt <- tau_dat %>% ggplot(aes(x = tau)) + geom_density() + theme_bw() + labs(x = "Posterior Tau",
                                                                               y = "Density") + 
    theme(text=element_text(size=7.5))
  return(plt)
}

plot_betas <- function(ghs_dat){
  beta_dat <- data.frame(selected = ifelse(ghs_dat$posterior_pe < 0.5, F, T),
                         j = ghs_dat$real_dat$signal$j, 
                         k = ghs_dat$real_dat$signal$k,
                         sig = ghs_dat$real_dat$signal$signal)
  plt <- beta_dat %>% ggplot(aes(x = as.factor(k), y = as.factor(j), color = as.factor(sig), fill = as.factor(sig), alpha = as.factor(selected))) + 
    geom_tile() + theme_bw() + labs(alpha = "Selection in Model",
                                    fill = "True Signal", 
                                    x = "Row Node",
                                    y = "Column Node") + guides(color = FALSE) + 
    theme(text=element_text(size=7.5))
  
}

plot_ghs <- function(ghs_dat, plot_graph = T){
  if(plot_graph)
    plot(ghs_dat$real_dat$full_ntwork)
  
  ((plot_lams(ghs_dat) +  plot_betas(ghs_dat)) / plot_gammas(ghs_dat)) + plot_tau(ghs_dat)
  
}

plot_comparison_lams_ghs <- function(ghs1, ghs2, names = NULL){
  if(!identical(ghs1$real_dat$signal, ghs2$real_dat$signal)){
    stop("Not the same network data/signal")
  }
  
  lam_dat1 <- data.frame(lam_mean = rowMeans(ghs1$posterior_draws_lam),
                         j = ghs1$real_dat$signal$j, 
                         k = ghs1$real_dat$signal$k,
                         sig = ghs1$real_dat$signal$signal
  )
  lam_dat2 <- data.frame(lam_mean = rowMeans(ghs2$posterior_draws_lam),
                         j = ghs2$real_dat$signal$j, 
                         k = ghs2$real_dat$signal$k,
                         sig = ghs2$real_dat$signal$signal
  )
  
  if(!is.null(names)){
    lam_dat_1 <- cbind.data.frame(lam_dat1, type = names[1] )
    lam_dat_2 <- cbind.data.frame(lam_dat2, type = names[2] )
  }
  
  connects_j <- lam_dat_1 %>% group_by(j) %>% summarise(j_cons = sum(sig))
  connects_k <- lam_dat_1 %>% group_by(k) %>% summarise(k_cons = sum(sig))
  
  lam_dat <- rbind.data.frame(lam_dat_1, lam_dat_2) %>% left_join(connects_j) %>% left_join(connects_k) %>% 
    mutate(cons = pmax(j_cons, k_cons))
  
  plt <- lam_dat %>% ggplot(aes(x = as.factor(k), y = as.factor(j), 
                                color = as.factor(sig), fill = as.factor(sig),
                                shape = as.factor(cons),
                                size = as.numeric(lam_mean))) + 
    geom_point() + theme_bw() + labs(size = "Mean Posterior Edge-Specific Shrinkage",
                                     color = "True Signal", 
                                     shape = "Connections",
                                     x = "Row Node",
                                     y = "Column Node") + guides(fill = FALSE) + 
    theme(text=element_text(size=7.5)) + facet_wrap(~type)
  
  return(plt)
}


