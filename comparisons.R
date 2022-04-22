seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#set.seed(505)
set.seed(seed)

#source the functions
source("GHSNetworkShrinkagePrior/generate_network_data.R")
source("GHSNetworkShrinkagePrior/posterior_computation.R")
source("GHSNetworkShrinkagePrior/compare_ghp.R")

burnins<-2500
draws<- 500
small_n <- 100
large_n <- 1000 
small_p <- 20
large_p <- 100

do_large_p <- F

#---
# Simulate Small P data =============================================================================================================
# --- 
small_n_small_p_data_no_hub <- generate_network_data(small_p,small_n,sigma2_beta = sqrt(25), prop_zero = 0.8)
large_n_small_p_data_no_hub <- generate_network_data(small_p,large_n,sigma2_beta = sqrt(25), prop_zero = 0.8)


small_num_nodes = sample(1:3,1)
small_nodes = sample(2:small_p - 3,small_num_nodes)
small_degrees = sample(1:10, length(small_nodes), replace = T)
small_n_small_p_data_hub <- generate_network_data(small_p,small_n,sigma2_beta = sqrt(15), prop_zero = 0.8,
                                                beta_hubs = T, hub_nodes = small_nodes, hub_degrees = small_degrees)
large_n_small_p_data_hub <- generate_network_data(small_p,large_n,sigma2_beta = sqrt(15), prop_zero = 0.8,
                                                 beta_hubs = T, hub_nodes = small_nodes, hub_degrees = small_degrees)



ghs_small_n_small_p_data_no_hub <- group_horseshoe_gibs(burnins,draws, dat = small_n_small_p_data_no_hub)
ghs_large_n_small_p_data_no_hub <- group_horseshoe_gibs(burnins,draws, dat = large_n_small_p_data_no_hub)
ghs_small_n_large_p_data_no_hub <- group_horseshoe_gibs(burnins,draws, dat = small_n_large_p_data_no_hub)
ghs_large_n_large_p_data_no_hub <- group_horseshoe_gibs(burnins,draws, dat = large_n_large_p_data_no_hub)

ghs_small_n_small_p_data_hub <- group_horseshoe_gibs(burnins,draws, dat = small_n_small_p_data_hub)
ghs_large_n_small_p_data_hub <- group_horseshoe_gibs(burnins,draws, dat = large_n_small_p_data_hub)
ghs_small_n_large_p_data_hub <- group_horseshoe_gibs(burnins,draws, dat = small_n_large_p_data_hub)
ghs_large_n_large_p_data_hub <- group_horseshoe_gibs(burnins,draws, dat = large_n_large_p_data_hub)
# ---
# Compare Small p data =========================================================================================================
# --- 
                        # NO HUBS # 
# SMALL N SMALL P ###
small_n_small_p_comparison_lm_nh <- compare_ghp(ghs_small_n_small_p_data_no_hub)
small_n_small_p_comparison_lasso_nh <- compare_ghp(ghs_small_n_small_p_data_no_hub, comparison = "lasso")
small_n_small_p_comparison_hp_nh <- compare_ghp(ghs_small_n_small_p_data_no_hub, comparison = "ungrouped")

small_n_small_p_no_hub <- list("LM" = small_n_small_p_comparison_lm_nh, 
                        "LASSO" = small_n_small_p_comparison_lasso_nh, 
                        "HP" = small_n_small_p_comparison_hp_nh, 
                        "GHSDat" = ghs_small_n_small_p_data_no_hub)

save(small_n_small_p_no_hub , file=paste0("/home/brennben/Robjects/", seed , "_small_n_small_p_no_hub.Rdata"))

                                      # LARGE N SMALL P 
large_n_small_p_comparison_lm_nh <- compare_ghp(ghs_large_n_small_p_data_no_hub)
large_n_small_p_comparison_lasso_nh <- compare_ghp(ghs_large_n_small_p_data_no_hub, comparison = "lasso")
large_n_small_p_comparison_hp_nh <- compare_ghp(ghs_large_n_small_p_data_no_hub, comparison = "ungrouped")

large_n_small_p_no_hub <- list("LM" = large_n_small_p_comparison_lm_nh, 
                               "LASSO" = large_n_small_p_comparison_lasso_nh, 
                               "HP" = large_n_small_p_comparison_hp_nh, 
                               "GHSDat" = ghs_large_n_small_p_data_no_hub)

save(large_n_small_p_no_hub , file=paste0("/home/brennben/Robjects/", seed , "_large_n_small_p_no_hub.Rdata"))

## HUBS ## 
# SMALL N SMALL P
small_n_small_p_comparison_lm_h <- compare_ghp(ghs_small_n_small_p_data_hub)
small_n_small_p_comparison_lasso_h <- compare_ghp(ghs_small_n_small_p_data_hub, comparison = "lasso")
small_n_small_p_comparison_hp_h <- compare_ghp(ghs_small_n_small_p_data_hub, comparison = "ungrouped")

small_n_small_p_hub <- list("LM" = small_n_small_p_comparison_lm_h, 
                            "LASSO" = small_n_small_p_comparison_lasso_h, 
                            "HP" = small_n_small_p_comparison_hp_h, 
                            "GHSDat" = ghs_small_n_small_p_data_hub)

save(small_n_small_p_hub , file=paste0("/home/brennben/Robjects/", seed , "_small_n_small_p_hub.Rdata"))

# LARGE N SMALL P 
large_n_small_p_comparison_lm_h <- compare_ghp(ghs_large_n_small_p_data_hub)
large_n_small_p_comparison_lasso_h <- compare_ghp(ghs_large_n_small_p_data_hub, comparison = "lasso")
large_n_small_p_comparison_hp_h <- compare_ghp(ghs_large_n_small_p_data_hub, comparison = "ungrouped")

large_n_small_p_hub <- list("LM" = large_n_small_p_comparison_lm_h, 
                            "LASSO" = large_n_small_p_comparison_lasso_h, 
                            "HP" = large_n_small_p_comparison_hp_h, 
                            "GHSDat" = ghs_large_n_small_p_data_hub)

save(large_n_small_p_hub , file=paste0("/home/brennben/Robjects/", seed , "_large_n_small_p_hub.Rdata"))

# ---
# ---
# ---
# ---
# ---
# ---
# ---
# ---
# Simulate Large P data =========================================================================================================
# --- 

if(do_large_p){
  
small_n_large_p_data_no_hub <- generate_network_data(large_p,small_n,sigma2_beta = sqrt(25), prop_zero = 0.9)
large_n_large_p_data_no_hub <- generate_network_data(large_p,large_n,sigma2_beta = sqrt(25), prop_zero = 0.9)

large_num_nodes = sample(1:15,1)
large_nodes = sample(25:large_p - 20,large_num_nodes)
large_degrees = sample(1:10, length(large_nodes), replace = T)
small_n_large_p_data_hub <- generate_network_data(large_p,small_n,sigma2_beta = sqrt(15), prop_zero = 0.8,
                                                  beta_hubs = T, hub_nodes = large_nodes, hub_degrees = large_degrees)
large_n_large_p_data_hub <- generate_network_data(large_p,large_n,sigma2_beta = sqrt(15), prop_zero = 0.8,
                                                  beta_hubs = T, hub_nodes = large_nodes, hub_degrees = large_degrees)


# ---
# Compare Large p Data =========================================================================================================
# --- 

# No Hubs 
# SMALL N LARGE P 
small_n_large_p_comparison_lm_nh <- compare_ghp(ghs_small_n_large_p_data_no_hub)
small_n_large_p_comparison_lasso_nh <- compare_ghp(ghs_small_n_large_p_data_no_hub, comparison = "lasso")
small_n_large_p_comparison_hp_nh <- compare_ghp(ghs_small_n_large_p_data_no_hub, comparison = "ungrouped")

small_n_large_p_no_hub <- list("LM" = small_n_large_p_comparison_lm_nh, 
                               "LASSO" = small_n_large_p_comparison_lasso_nh, 
                               "HP" = small_n_large_p_comparison_hp_nh, 
                               "GHSDat" = ghs_small_n_large_p_data_no_hub)
save(small_n_large_p_no_hub , file=paste0("/home/brennben/Robjects/", seed , "_small_n_large_p_no_hub.Rdata"))

# LARGE N LARGE P 
large_n_large_p_comparison_lm_nh <- compare_ghp(ghs_large_n_large_p_data_no_hub)
large_n_large_p_comparison_lasso_nh <- compare_ghp(ghs_large_n_large_p_data_no_hub, comparison = "lasso")
large_n_large_p_comparison_hp_nh <- compare_ghp(ghs_large_n_large_p_data_no_hub, comparison = "ungrouped")

large_n_large_p_no_hub <- list("LM" = large_n_large_p_comparison_lm_nh, 
                               "LASSO" = large_n_large_p_comparison_lasso_nh, 
                               "HP" = large_n_large_p_comparison_hp_nh, 
                               "GHSDat" = ghs_large_n_large_p_data_no_hub)

save(large_n_large_p_no_hub , file=paste0("/home/brennben/Robjects/", seed , "_large_n_large_p_no_hub.Rdata"))

# Hubs 

# SMALL N LARGE P 
small_n_large_p_comparison_lm_h <- compare_ghp(ghs_small_n_large_p_data_hub)
small_n_large_p_comparison_lasso_h <- compare_ghp(ghs_small_n_large_p_data_hub, comparison = "lasso")
small_n_large_p_comparison_hp_h <- compare_ghp(ghs_small_n_large_p_data_hub, comparison = "ungrouped")

small_n_large_p_hub <- list("LM" = small_n_large_p_comparison_lm_h, 
                               "LASSO" = small_n_large_p_comparison_lasso_h, 
                               "HP" = small_n_large_p_comparison_hp_h, 
                               "GHSDat" = ghs_small_n_large_p_data_hub)
save(small_n_large_p_hub , file=paste0("/home/brennben/Robjects/", seed , "_small_n_large_p_hub.Rdata"))

# LARGE N LARGE P 
large_n_large_p_comparison_lm_h <- compare_ghp(ghs_large_n_large_p_data_hub)
large_n_large_p_comparison_lasso_h <- compare_ghp(ghs_large_n_large_p_data_hub, comparison = "lasso")
large_n_large_p_comparison_hp_h <- compare_ghp(ghs_large_n_large_p_data_hub, comparison = "ungrouped")

large_n_large_p_hub <- list("LM" = large_n_large_p_comparison_lm_h, 
                               "LASSO" = large_n_large_p_comparison_lasso_h, 
                               "HP" = large_n_large_p_comparison_hp_h, 
                               "GHSDat" = ghs_large_n_large_p_data_hub)
save(large_n_large_p_hub , file=paste0("/home/brennben/Robjects/", seed , "_large_n_large_p_hub.Rdata"))

} 


