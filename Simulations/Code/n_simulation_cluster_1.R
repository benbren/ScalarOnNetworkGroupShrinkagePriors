seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))



#source the functions
source("GHSNetworkShrinkagePrior/generate_network_data.R")
source("GHSNetworkShrinkagePrior/posterior_computation.R")
source("GHSNetworkShrinkagePrior/compare_ghp.R")



sims_MSE <- NULL 
sims_signalMSE <- NULL 
sims_zeroMSE <- NULL 
sims_accuracy <- NULL 
sims_pred <- NULL 
sims_sens <- NULL 
sims_spec <- NULL

for (ss in c(100, 500, 1000, 2500,5000)) {
  
  set.seed(seed)
  
  dat <- generate_network_data(p = 50, n = ss, r2 = 0.25, groups = T, num_groups = 5, 
                               group_sizes = rep(10,5), num_signals = c(0,0,45,rep(0,2)), 
                               beta_mean = 5, beta_sd = 1)
  ghs <- group_horseshoe_gibs_groups(1250,250, dat = dat, init_beta = 10)
  save(ghs , file=paste0("/home/brennben/Robjects/", seed , "ghs", ss, ".RData"))
  
  unghp <- compare_ghp(ghs, comparison = "ungrouped")
  
  sims_MSE <- cbind(sims_MSE, c(unghp$ghp$MSE,unghp$comparison$MSE))
  sims_signalMSE <- cbind(sims_signalMSE, c(unghp$ghp$SignalMSE,unghp$comparison$SignalMSE))
  sims_zeroMSE <- cbind(sims_zeroMSE, c(unghp$ghp$ZeroMSE,unghp$comparison$ZeroMSE))
  sims_sens <- cbind(sims_sens, c(unghp$ghp$Sensitivity,unghp$comparison$Sensitivity))
  sims_spec <- cbind(sims_spec, c(unghp$ghp$Specificity, unghp$comparison$Specificity))
  sims_accuracy <- cbind(sims_accuracy, c(unghp$ghp$Accuracy,unghp$comparison$Accuracy))
  sims_pred <- cbind(sims_pred, c(unghp$ghp$PredictionError,unghp$comparison$PredictionError))
  
}


sim_dat <- rbind(sims_MSE, sims_signalMSE, sims_zeroMSE,sims_accuracy, sims_sens,sims_spec,sims_pred)
colnames(sim_dat) <- as.character(c(100, 500, 1000, 2500,5000))
sim_dat <- as.data.frame(sim_dat)
sim_dat$Metric <- c(rep("MSE", 2),rep("SignalMSE", 2), rep("ZeroMSE", 2),rep("Accuracy",2), rep("Sensitivity", 2), rep("Specificity", 2), rep("Prediction Error", 2))
sim_dat$Model <- rep(c("GHP", "UnGHP"), 7)


sim_dat_long <- sim_dat %>% pivot_longer(-c(Metric, Model), names_to = "N")

save(sim_dat_long , file=paste0("/home/brennben/Robjects/", seed , "_n_simulation.Rdata"))