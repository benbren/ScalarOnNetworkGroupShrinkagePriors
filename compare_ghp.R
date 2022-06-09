bayesian_variable_selection <- function(ghs_model, alpha, delta = NULL){
  
  s <- ghs_model$real_dat$signal$signal == 1
  if(is.null(delta)){
    delta <- 0.5*min(ghs_model$real_dat$signal$beta[s])
  }
  
  choices <- rep(0,choose(ghs_model$real_dat$p,2))
  ps <- rowMeans(abs(ghs_model$posterior_draws_beta) > delta)
  ps < ifelse(ps == 1, 1 - ((2*ghs_model$posterior_draws)^-1), ps)
  ps_sorted <- sort(ps, decreasing = T)
  phi <- 0 
  for (i in 1:length(ps_sorted)){
    if(mean(1-ps_sorted[1:i]) <= alpha){ 
      next
    } else{
      phi = ps_sorted[i-1]
      break
    }
  }
  choices_idx <- which(ps >= phi)
  choices[choices_idx] <- 1 
  
  return(choices)
}

inclusion_weight_selection <- function(ghs_model){
  choices <- rep(0,choose(ghs_model$real_dat$p,2))
  tau2 <- mean(ghs_model$posterior_draws_tau)
  connection_specific <- rowMeans(ghs$posterior_draws_lam)*rowMeans(ghs$posterior_draws_gam)
  w <- i / (1 + connection_specific*tau2)
  d <- which(w > 0.5)
  choices[d] <- 1 
  return(choices)
}

summarize_ghs_metrics <- function(ghs_model, fdr = 0.05){
  
  n <- ghs_model$real_dat$n
  p <- ghs_model$real_dat$p
  q <- choose(p,2)
  real_betas <- ghs_model$real_dat$signal$beta
  signals <- ghs_model$real_dat$signal$signal == 1
  y_true <- ghs_model$real_dat$outcomes
  
  mse_ghp <- sum((real_betas - ghs_model$posterior_pe)^2) / q
  mse_ghp_s <- sum((real_betas[signals] - ghs_model$posterior_pe[signals])^2) / sum(signals)
  mse_ghp_z <- sum((real_betas[!signals] - ghs_model$posterior_pe[!signals])^2) / sum(!signals)
  ghp_choices <- bayesian_variable_selection(ghs_model, alpha = fdr)
  ghp_accuracy = sum(ghp_choices == ghs_model$real_dat$signal$signal)/ q
  ghp_sens = sum(ghp_choices[signals] == 1)/ sum(signals)
  ghp_spec = sum(ghp_choices[!signals] == 0)/ sum(!signals)

  ghp_preds <- ghs_model$predictions
  pred_error_ghp <- sum((y_true - ghp_preds)^2) / n
  
  return(list(MSE = mse_ghp, MSESignal = mse_ghp_s, MSEZero = mse_ghp_z, 
              ShrinkageAccuracy = ghp_accuracy, Sensitivity = ghp_sens, Specificity = ghp_spec,
              PredictionError = pred_error_ghp))
} 


compare_ghp <- function(ghs_model,comparison = "lm", nodes = F, delta = NULL, alpha = 0.1){
  
  n <- ghs_model$real_dat$n
  p <- ghs_model$real_dat$p
  q <- choose(p,2)

  X <- t(matrix(unlist(ghs_model$real_dat$ntworkslt), choose(p,2), n))
  y_true <- ghs_model$real_dat$outcomes
  
  mod_dat <- as.data.frame(cbind(y_true, X))
  names(mod_dat) <- c("y", paste0("X",1:dim(X)[2]))
  real_betas <- ghs_model$real_dat$signal$beta
  signals <- ghs_model$real_dat$signal$signal == 1
  if(is.null(delta)){
    delta <- 0.5*min(abs(ghs_model$real_dat$signal$beta[signals]))
  }

  if(comparison == "lm"){
    
    print("Linear Model")
    
    comparison_mod <- lm(y ~ ., dat = mod_dat)
    comparison_coefficients <- comparison_mod$coefficients[-1]
    mse_comparison <- sum((real_betas - comparison_coefficients)^2) / q
    mse_comparison_s <- sum((real_betas[signals] - comparison_coefficients[signals])^2) / sum(signals)
    mse_comparison_z <- sum((real_betas[!signals] - comparison_coefficients[!signals])^2) / sum(!signals)
    
    comparison_choices <- rep(0,choose(p,2))
    comparison_choices_idx <- which(summary(comparison_mod)$coefficients[-1,4] < 0.05)
    comparison_choices[comparison_choices_idx] <- 1
    comparison_accuracy <- sum(comparison_choices == ghs_model$real_dat$signal$signal)/ q
    comparison_sens = sum(comparison_choices[signals] == 1)/ sum(signals)
    comparison_spec = sum(comparison_choices[!signals] == 0)/ sum(!signals)
    
    comparison_preds <- predict(comparison_mod)
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n

    
  } else if(comparison == "lasso"){
    print("LASSO")
    
    library(glmnet)
    comparison_mod <- glmnet::glmnet(X, y_true, lambda = 0.1,
                             family = "gaussian", 
                             alpha = 1)
    
    comparison_coefficients <- comparison_mod$beta
    mse_comparison <- sum((real_betas - comparison_coefficients)^2) / q
    mse_comparison_s <- sum((real_betas[signals] - comparison_coefficients[signals])^2) / sum(signals)
    mse_comparison_z <- sum((real_betas[!signals] - comparison_coefficients[!signals])^2) / sum(!signals)
    
    comparison_choices <- rep(0,q)
    comparison_choices_idx <- comparison_mod$beta@i
    comparison_choices[comparison_choices_idx] <- 1
    comparison_accuracy = sum(comparison_choices == ghs_model$real_dat$signal$signal)/ q
    comparison_sens = sum(comparison_choices[signals] == 1)/ sum(signals)
    comparison_spec = sum(comparison_choices[!signals] == 0)/ sum(!signals)
    comparison_preds <- X%*%comparison_coefficients
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n
  
  } else if(comparison == "ungrouped"){
    
    if(nodes){
      
      comparison_mod <- group_horseshoe_gibs_nodes(burn_ins = ghs_model$burn_ins,
                                             posterior_draws = ghs_model$posterior_draws,
                                             dat = ghs_model$real_dat, ungroup = T, generate = F)
    } else {
      comparison_mod <- group_horseshoe_gibs_groups(burn_ins = ghs_model$burn_ins,
                                             posterior_draws = ghs_model$posterior_draws,
                                             dat = ghs_model$real_dat, ungroup = T)
    } 
    
    comparison_coefficients <- comparison_mod$posterior_pe

    mse_comparison <- sum((real_betas - comparison_mod$posterior_pe)^2) / q
    mse_comparison_s <- sum((real_betas[signals] - comparison_coefficients[signals])^2) / sum(signals)
    mse_comparison_z <- sum((real_betas[!signals] - comparison_coefficients[!signals])^2) / sum(!signals)
    comparison_choices <- inclusion_weight_selection(comparison_mod, alpha)
    comparison_accuracy = sum(comparison_choices == comparison_mod$real_dat$signal$signal)/ q
    comparison_sens = sum(comparison_choices[signals] == 1)/ sum(signals)
    comparison_spec = sum(comparison_choices[!signals] == 0)/ sum(!signals)
    comparison_preds <- comparison_mod$predictions
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n
  } 
  
 
  print("GHP")
  
  ghp_coefficients <- ghs_model$posterior_pe
  
  mse_ghp <- sum((real_betas - ghp_coefficients)^2) / q
  mse_ghp_s <- sum((real_betas[signals] - ghp_coefficients[signals])^2) / sum(signals)
  mse_ghp_z <- sum((real_betas[!signals] - ghp_coefficients[!signals])^2) / sum(!signals)
  ghp_choices <- inclusion_weight_selection(ghs_model,alpha)
  ghp_accuracy = sum(ghp_choices == ghs_model$real_dat$signal$signal)/ q
  ghp_sens = sum(ghp_choices[signals] == 1)/ sum(signals)
  ghp_spec = sum(ghp_choices[!signals] == 0)/ sum(!signals)
  ghp_preds <- ghs_model$predictions
  pred_error_ghp <- sum((y_true - ghp_preds)^2) / n
  
  #return(return_dat)
  
  return(list(comparison = list(MSE = mse_comparison, SignalMSE = mse_comparison_s, ZeroMSE = mse_comparison_z,
                                Accuracy = comparison_accuracy, PredictionError = pred_error_comparison,
                                Sensitivity = comparison_sens, Specificity = comparison_spec),
              ghp =  list(MSE = mse_ghp, SignalMSE = mse_ghp_s, ZeroMSE = mse_ghp_z,
                          Accuracy = ghp_accuracy, PredictionError = pred_error_ghp,
                          Sensitivity = ghp_sens, Specificity = ghp_spec)))
  
}






