################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from linear model G-computation
# 18/06/2025
################################################################################


### Function to do G-computation with linear model
#-------------------------------------------------
run_lm_gcomp_boot <- function(obs, exposure_shift = 1, n_boot = 200, alpha = 0.05, sim_id = 1, psi) {
  
  
  outcome_name <- names(obs)[1]  # "Y"
  covariate_start <- which(names(obs) == "age")[1]
  exposures <- names(obs)[2:(covariate_start - 1)]
  covariates <- names(obs)[covariate_start:length(obs)]
  
  
  ### Function to compute ATE for a given data set
  #-----------------------------------------------
  compute_ate <- function(data) {
    formula_rf <- as.formula(paste(outcome_name, "~", paste(c(exposures, covariates), collapse = " + ")))
    
    LM_model <- lm(formula_rf, data = data)
    
    Y_hat_observed <- predict(LM_model, newdata = data)
    
    data_shifted <- data
    data_shifted[, exposures] <- data_shifted[, exposures] + exposure_shift
    
    Y_hat_shifted <- predict(LM_model, newdata = data_shifted)
    
    mean(Y_hat_shifted - Y_hat_observed)
  }
  
  
  ### Estimate ATE on original sample
  #----------------------------------
  ate_point_estimate <- compute_ate(obs)
  
  
  ### Bootstrap
  #------------
  ate_boot <- numeric(n_boot)
  for (b in 1:n_boot) {
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_sample <- obs[boot_indices, ]
    ate_boot[b] <- compute_ate(boot_sample)
  }
  
  
  ### Calculate bootstrap statistics
  #---------------------------------
  ci_lower <- quantile(ate_boot, alpha / 2)
  ci_upper <- quantile(ate_boot, 1 - alpha / 2)
  empirical_power <- (ci_lower > 0 | ci_upper < 0)
  CI_coverage <- (ci_lower <= psi[1]) & (ci_upper >= psi[1])
  
  
  ### Format output
  #----------------
  coef_rf <- data.frame(
    Estimate = ate_point_estimate,
    Sim = as.factor(sim_id),
    parameter = "Sum betas",
    Method = paste0("Linear model g-computation (bs=", n_boot, ")"),
    Size = n,
    Power = empirical_power,
    CI_coverage = CI_coverage,
    CI_width = ci_upper-ci_lower,
    PIP = NA
  )
  
  return(coef_rf)
}

