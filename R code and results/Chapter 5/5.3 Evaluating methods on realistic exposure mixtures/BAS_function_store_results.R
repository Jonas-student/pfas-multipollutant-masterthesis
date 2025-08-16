################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from Bayesian additive sampling
# 17/06/2025
################################################################################


### Load libraries 
#-----------------
library(BAS)


### Function to store BAS results
#--------------------------------
run_BAS <- function(obs, alpha = 0.05, sim_id = 1, psi){
    
    
  ### Fit BAS model
  #----------------
  bas_fit <- bas.lm(
    Y ~ lbpfos + pfba + pfda + lbpfhxs + pfna + lbpfoa + pfos, 
    include.always = as.formula(~age + gender + Educ2 + Educ3),
    data = obs,
    method="MCMC+BAS", 
    n.models=2^10, 
    MCMC.iterations=50000, 
    thin=10, 
    prior="JZS",
    modelprior = uniform()
  )
    
    
  # Extract model-averaged coefficients
  coef_summary <- coef(bas_fit, estimator = "BMA")
    
    
  # Coefficients
  chemical_effects <- coef_summary$postmean[2:8]
  exposure_ses <- coef_summary$postsd[2:8]
    
    
  # Save effect
  joint_effects <- sum(chemical_effects)
  
  
  ### Joint effect summary
  #-----------------------
  ci_width <- NA
  empirical_power <- NA
  CI_coverage <- NA
  
  
  ### Define all parameter names and corresponding estimates
  #---------------------------------------------------------
  parameters <- c("Joint_Increase_All_Exposures", "lbpfos", 
                  "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
  
  
  ### Placeholder vectors
  #----------------------
  chemical_estimates <- numeric(7)
  ci_lowers <- numeric(7)
  ci_uppers <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  
  
  ### Calculate statistics for each chemical
  #-----------------------------------------
  for(j in seq(1,7,1)){
    
    chemical_estimates[j] <- chemical_effects[j]
    
    z <- qnorm(1 - alpha/2)
    ci <- c(
      chemical_effects[j] - z * exposure_ses[j],
      chemical_effects[j] + z * exposure_ses[j]
    )
    
    ci_widths[j] <- ci[2] - ci[1]
    empirical_powers[j] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[j] <- (ci[1] <= psi[1] * psi[j+1]) & (ci[2] >= psi[1] * psi[j+1])
  }
  
  
  ### Add joint effect to the chemical list
  #----------------------------------------
  estimates <- c(joint_effects, chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Store and return as a dataframe
  #----------------------------------
  coef_bayes <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = paste0("Bayesian model averaging"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all, 
    PIP = c(NA , coef_summary$probne0[-1])
  )
  
  return(coef_bayes)
}
