compare_ghp <- function(horseshoe_data,comparison = "lm"){
  
  n <- horseshoe_data$real_dat$n
  p <- horseshoe_data$real_dat$p
  q <- choose(p,2)

  X <- t(matrix(unlist(horseshoe_data$real_dat$ntworkslt), choose(p,2), n))
  y_true <- horseshoe_data$real_dat$outcomes
  
  mod_dat <- as.data.frame(cbind(y_true, X))
  names(mod_dat) <- c("y", paste0("X",1:dim(X)[2]))
  real_betas <- horseshoe_data$real_dat$signal$beta

  if(comparison == "lm"){
    
    print("Linear Model")
    
    comparison_mod <- lm(y ~ ., dat = mod_dat)
    comparison_coefficients <- comparison_mod$coefficients[-1]
    mse_comparison <- sum((real_betas - comparison_coefficients)^2) / q
    
    comparison_choices <- rep(0,choose(p,2))
    comparison_choices_idx <- which(summary(comparison_mod)$coefficients[-1,4] < 0.05)
    comparison_choices[comparison_choices_idx] <- 1
    comparison_accuracy = sum(comparison_choices == horseshoe_data$real_dat$signal$signal)/ q
    
    comparison_preds <- predict(comparison_mod)
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n

    
  } else if(comparison == "lasso"){
    print("LASSO")
    
    library(glmnet)
    comparison_mod <- glmnet(X, y_true, lambda = 0.1,
                             family = "gaussian", 
                             alpha = 1)
    
    comparison_coefficients <- comparison_mod$beta
    mse_comparison <- sum((real_betas - comparison_coefficients)^2) / q
    
    comparison_choices <- rep(0,q)
    comparison_choices_idx <- comparison_mod$beta@i
    comparison_choices[comparison_choices_idx] <- 1
    comparison_accuracy = sum(comparison_choices == horseshoe_data$real_dat$signal$signal)/ q
    comparison_preds <- X%*%comparison_coefficients
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n
  
  } else if(comparison == "ungrouped"){
    print("HP")
    
    comparison_mod <- group_horseshoe_gibs(burn_ins = horseshoe_data$burn_ins,
                                           posterior_draws = horseshoe_data$posterior_draws,
                                           dat = horseshoe_data$real_dat, ungroup = T, generate = F)
    
    mse_comparison <- sum((real_betas - comparison_mod$posterior_pe)^2) / q
    
    comparison_choices <- rep(0,q)
    comparison_choices_idx <- which(comparison_mod$posterior_pe - 1.96*comparison_mod$posterior_se > 0 | comparison_mod$posterior_pe + 1.96*comparison_mod$posterior_se < 0)
    comparison_choices[comparison_choices_idx] <- 1 
    comparison_accuracy = sum(comparison_choices == comparison_mod$real_dat$signal$signal)/ q
    
    comparison_preds <- comparison_mod$predictions
    pred_error_comparison <- sum((y_true - comparison_preds)^2) / n
  } 
  
 
  print("GHP")
  
  mse_ghp <- sum((real_betas - horseshoe_data$posterior_pe)^2) / q
  
  ghp_choices <- rep(0,q)
  ghp_choices_idx <- which(horseshoe_data$posterior_pe - 1.96*horseshoe_data$posterior_se > 0 | horseshoe_data$posterior_pe + 1.96*horseshoe_data$posterior_se < 0)
  ghp_choices[ghp_choices_idx] <- 1 
  ghp_accuracy = sum(ghp_choices == horseshoe_data$real_dat$signal$signal)/ q
  
  ghp_preds <- horseshoe_data$predictions
  pred_error_ghp <- sum((y_true - ghp_preds)^2) / n
  
  #return(return_dat)
  
  return(list(comparison = list(MSE = mse_comparison, Accuracy = comparison_accuracy, PredictionError = pred_error_comparison),
              ghp =  list(MSE = mse_ghp, Accuracy = ghp_accuracy, PredictionError = pred_error_ghp)))
  
}






