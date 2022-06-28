dat <- generate_network_data(p = 30, n = 500, r2 = noise, groups = T, num_groups = 5, 
                             group_sizes = rep(6,5), num_signals = c(0,0,15,0,0), 
                             beta_mean = 5, beta_sd = 1) 

gh <- rstan::stan(file  = "group_horseshoe.stan", data = dat$stan_dat, chains = 4L, iter = 4000L, warmup = 3000L)

gh_fit <- rstan::stan(fit = gh, data = stan_dat,  chains = 2, iter = 2000, refresh = 0)
   