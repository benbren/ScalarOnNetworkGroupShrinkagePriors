summarize_group_boosting <- function(dat){
  test = GroupBoosting::group_boosting_fit(t(matrix(unlist(dat$ntworkslt), choose(dat$p,2), dat$n)), 
                                           dat$outcomes, 
                                           dat$signal$group)
  s = dat$signal$signal == 1
  mse = mean((test$beta - dat$signal$beta)^2)
  zeromse = mean((test$beta[!s] - dat$signal$beta[!s])^2)
  signalmse = mean((test$beta[s] - dat$signal$beta[s])^2)
  group_boosting_choices = ifelse(test$beta == 0, 0, 1)
  accuracy = mean(group_boosting_choices == dat$signal$signal)
  sensitivity = mean(group_boosting_choices[s] == dat$signal$signal[s])
  specificity = mean(group_boosting_choices[!s] == dat$signal$signal[!s])
  return(data.frame(Model = "GroupBoosting", 
                    Metric = c("MSE", "ZeroMSE", "SignalMSE", "Accuracy", "Sensitivity", "Specificity"),
                    value = c(mse, zeromse, signalmse, accuracy, sensitivity, specificity)))
}

gb_add_on <- NULL
for (ss in c(100, 500,1000,2500)){
  print(ss)
  set.seed(505)
  dat <- generate_network_data(p = 50, n = ss, r2 = 0.25, groups = T, num_groups = 5, 
                               group_sizes = rep(10,5), num_signals = c(0,0,45,rep(0,2)), 
                               beta_mean = 5, beta_sd = 1)
  sim_dat <- cbind(summarize_group_boosting(dat), ss)
  gb_add_on <- rbind(gb_add_on, sim_dat)
  
}


models <- c("GHP","UnGHP", "GroupBoosting")
metrics <- c("MSE", "ZeroMSE", "SignalMSE")
metrics1 <- c("Accuracy", "Sensitivity", "Specificity")

sim_dat_n %>% filter(Model %in% models & Metric %in% metrics1)  %>% filter(N!= 5000) %>%  mutate(N = as.numeric(N))%>%  
  ggplot(aes(x = N, y = value, group = Model, 
                            color = Model)) + geom_point() + 
  geom_line() + facet_grid(Metric~., scales = "free_y") + theme_light()

sim_dat_n %>% filter(Metric == "ZeroMSE") %>% pivot_wider(id_cols = N, names_from = Model)

