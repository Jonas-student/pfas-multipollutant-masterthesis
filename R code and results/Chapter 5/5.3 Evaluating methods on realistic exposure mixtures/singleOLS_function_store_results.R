################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from single OLS regression model
# 24/06/2025
################################################################################


run_singlelinear_model <- function(obs, alpha = 0.05, sim_id = 1, psi) {

  
  ### Placeholders for chemical results
  #------------------------------------
  chemical_estimates <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  exposure_names<-c("lbpfos", 
    "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    
    
    ### Fit single-pollutant linear model
    #------------------------------------
    formula <- as.formula(paste0("Y ~ age + gender + Educ2 + Educ3 + ", exposure_names[i]))
    linear_model <- lm(formula, data = obs)
    
    
    ### Select the estimate
    #----------------------
    chemical_estimates[i] <- coef(linear_model)[exposure_names[i]]
    
    
    ### Select the CI parameters
    #---------------------------
    ci <- confint(linear_model)[exposure_names[i],]
    ci_widths[i] <- ci[2] - ci[1]
    empirical_powers[i] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[i] <- (ci[1] <= psi[1]*psi[i + 1]) & (ci[2] >= psi[1]*psi[i + 1])
  }
  
  ### Prepare output
  #-----------------
  coef_lm <- data.frame(
    Estimate = chemical_estimates,
    Sim = as.factor(sim_id),
    Size = nrow(obs),
    parameter = exposure_names,
    Method = "Single pollutant linear model",
    Power = empirical_powers,
    CI_coverage = CI_coverages,
    CI_width = ci_widths,
    PIP = rep(NA,7)
  )
  
  return(coef_lm)
}

