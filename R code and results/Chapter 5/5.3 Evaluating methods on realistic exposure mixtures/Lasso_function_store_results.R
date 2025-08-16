################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from Lasso regression
# 05/07/2025
################################################################################


### Load libraries 
#-----------------
library(glmnet)


### Function to store BAS results
#--------------------------------
run_Lasso <- function(obs, n_boot = 200, alpha = 0.05, sim_id = 1, psi){
  
  
  ### Prepare data
  #---------------
  X <- model.matrix(Y ~ age + gender + Educ2 + Educ3 + lbpfos + pfba + pfda + lbpfhxs + pfna + lbpfoa + pfos, data = obs)[, -1]
  Y <- obs$Y
  penalty_factors <- c(rep(0, 4), rep(1, 7))
  
  
  ### Store bootstrap results
  #--------------------------
  joint_effects <- numeric(n_boot)
  chemical_effects <- matrix(nrow = n_boot, ncol = 7)
  chemical_selected <- matrix(FALSE, nrow = n_boot, ncol = 7)
  
  
  ### Run bootstrap sampling
  #-------------------------
  for(i in 1:n_boot){
    
    
    # Bootstrap sample
    boot_indices <- sample(1:nrow(obs), replace = TRUE)
    X_boot <- X[boot_indices, ]
    Y_boot <- Y[boot_indices]
    
    
    # Fit Ridge regression using cross-validated lambda
    cv_fit <- cv.glmnet(X_boot, Y_boot, alpha = 1, penalty.factor = penalty_factors,
                        type.measure = "mse", nfolds = 5)
    lasso_fit <- glmnet(X_boot, Y_boot, alpha = 1, penalty.factor = penalty_factors,
                        lambda = cv_fit$lambda.1se)
    
    
    # Extract coefficients
    coef_lasso <- as.vector(coef(lasso_fit))
    
    
    # Extract exposure estimates (positions 6 to 12: lbpfos to pfos)
    exposure_estimates <- coef_lasso[6:12]
    
    
    # Save joint and individual effects
    joint_effects[i] <- sum(exposure_estimates)
    chemical_effects[i, ] <- exposure_estimates
    
    
    # Mark selected if coef != 0
    chemical_selected[i, ] <- exposure_estimates != 0
  }
  
  
  ### Joint effect summary
  #-----------------------
  ci <- quantile(joint_effects, probs = c(0.025, 0.975))
  ci_lower <- ci[1]
  ci_upper <- ci[2]
  ci_width <- ci_upper - ci_lower
  empirical_power <- (ci_lower > 0 | ci_upper < 0)
  CI_coverage <- (ci_lower <= psi[1]) & (ci_upper >= psi[1])
  
  
  ### Parameter names
  #------------------
  chemicals <- c("lbpfos", "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
  
  
  ### Placeholders for chemical results
  #------------------------------------
  chemical_estimates <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  selection_freqs <- numeric(7)   # <-- PIP alternative
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    effect_samples <- chemical_effects[, i]
    
    chemical_estimates[i] <- mean(effect_samples)
    
    ci <- quantile(effect_samples, probs = c(0.025, 0.975))
    ci_widths[i] <- ci[2] - ci[1]
    empirical_powers[i] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[i] <- (ci[1] <= psi[1]*psi[i + 1]) & (ci[2] >= psi[1]*psi[i + 1])
    
    # Calculate proportion of times selected (nonzero coef)
    selection_freqs[i] <- mean(chemical_selected[, i])
  }
  
  
  ### Store all results
  #--------------------
  parameters <- c("Joint_Increase_All_Exposures", chemicals)
  estimates <- c(mean(joint_effects), chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Store results together
  #-------------------------
  coef_lasso <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = paste0("Lasso regression (bs=", n_boot, ")"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = c(NA, selection_freqs)
  )
  
  
  return(coef_lasso)
}
