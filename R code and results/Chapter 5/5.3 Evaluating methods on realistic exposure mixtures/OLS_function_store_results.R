################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from OLS regression model
# 17/06/2025
################################################################################


run_linear_model <- function(obs, alpha = 0.05, sim_id = 1, psi) {
  
  
  ### Fit linear model
  #-------------------
  linear_model <- lm(Y ~ age + gender + Educ2 + Educ3 + lbpfos + pfba + pfda + lbpfhxs + pfna + lbpfoa + pfos, data = obs)
  
  
  ### Covariance matrix
  #--------------------
  cov_matrix <- vcov(linear_model)
  
  
  ### Exposure names
  #-----------------
  exposure_names <- c("lbpfos", "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
  
  
  ### Extract coefficients
  #-----------------------
  coefs <- coef(linear_model)[exposure_names]
  
  
  ### Sum of coefficients (joint effect)
  #-------------------------------------
  sum_coef <- sum(coefs)
  
  
  ### Variance and standard error of the sum
  #-----------------------------------------
  var_sum <- sum(diag(cov_matrix[exposure_names, exposure_names])) + 
    2 * sum(cov_matrix[exposure_names, exposure_names][lower.tri(cov_matrix[exposure_names, exposure_names])])
  se_sum <- sqrt(var_sum)
  
  
  ### Calculate confidence interval
  #--------------------------------
  ci_lower <- sum_coef - qnorm(1 - alpha/2) * se_sum
  ci_upper <- sum_coef + qnorm(1 - alpha/2) * se_sum
  ci_width <- ci_upper - ci_lower
  
  
  ### Calculate power and CI coverage
  #----------------------------------
  z_value <- sum_coef / se_sum
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  empirical_power <- p_value < alpha
  CI_coverage <- (ci_lower <= psi[1]) & (ci_upper >= psi[1])
  
  
  ### Placeholders for chemical results
  #------------------------------------
  chemical_estimates <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    
    chemical_estimates[i] <- coef(linear_model)[exposure_names[i]]
    
    ci <- confint(linear_model)[exposure_names[i],]
    ci_widths[i] <- ci[2] - ci[1]
    empirical_powers[i] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[i] <- (ci[1] <= psi[1]*psi[i + 1]) & (ci[2] >= psi[1]*psi[i + 1])
  }
  
  
  ### Store all results
  #--------------------
  parameters <- c("Joint_Increase_All_Exposures", exposure_names)
  estimates <- c(sum_coef, chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Prepare output
  #-----------------
  coef_lm <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    Size = nrow(obs),
    parameter = parameters,
    Method = "Multiple pollutant linear model",
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = rep(NA,8)
  )
  
  return(coef_lm)
}

