seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


set.seed(seed)

#source the functions
source("GHSNetworkShrinkagePrior/generate_network_data.R")
source("GHSNetworkShrinkagePrior/posterior_computation.R")
source("GHSNetworkShrinkagePrior/compare_ghp.R")



sims_MSE <- NULL 
sims_accuracy <- NULL 
sims_pred <- NULL 
sims_sens <- NULL 
sims_spec <- NULL

set.seed(15)

for (noise in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  
  dat <- generate_network_data(p = 75, n = 500, r2 = noise, groups = T, num_groups = 7, 
                               group_sizes = c(14,5,11,12,8,10,15), num_signals = c(0,0,12,0,0,0,0,5,0), 
                               beta_mean = 3, beta_sd = 1)
  ghs <- group_horseshoe_gibs_groups(1250,250, dat = dat)
  lin <- compare_ghp(ghs)
  lasso <-compare_ghp(ghs, comparison = "lasso")
  unghp <- compare_ghp(ghs, comparison = "ungrouped")
  
  sims_MSE <- cbind(sims_MSE, c(lin$ghp$MSE, lin$comparison$MSE, lasso$comparison$MSE, unghp$comparison$MSE))
  sims_sens <- cbind(sims_sens, c(lin$ghp$Sensitivity, lin$comparison$Sensitivity, lasso$comparison$Sensitivity, unghp$comparison$Sensitivity))
  sims_spec <- cbind(sims_spec, c(lin$ghp$Specificity, lin$comparison$Specificity, lasso$comparison$Specificity, unghp$comparison$Specificity))
  sims_accuracy <- cbind(sims_accuracy, c(lin$ghp$Accuracy, lin$comparison$Accuracy, lasso$comparison$Accuracy, unghp$comparison$Accuracy))
  sims_pred <- cbind(sims_pred, c(lin$ghp$PredictionError, lin$comparison$PredictionError, lasso$comparison$PredictionError, unghp$comparison$PredictionError))
  
}


sim_dat <- rbind(sims_MSE, sims_accuracy, sims_sens,sims_spec,sims_pred)
colnames(sim_dat) <- as.character(c(0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8, 0.9))
sim_dat <- as.data.frame(sim_dat)
sim_dat$Metric <- c(rep("MSE", 4), rep("Accuracy",4), rep("Sensitivity", 4), rep("Specificity", 4), rep("Prediction Error", 4))
sim_dat$Model <- rep(c("GHP", "LM", "LASSO", "UnGHP"), 5)


sim_dat_long <- sim_dat %>% pivot_longer(-c(Metric, Model), names_to = "StN Ratio")

save(sim_dat_long , file=paste0("/home/brennben/Robjects/", seed , "_stn_simulation.Rdata"))