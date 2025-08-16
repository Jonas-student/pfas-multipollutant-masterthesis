################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from Bayesian kernel machine regression
# 17/06/2025
################################################################################


### Load libraries 
#-----------------
library(bkmr)


### Function to store BAS results
#--------------------------------
run_bkmr <- function(obs, sim_id = 1, psi){
  
  
  ### Prepare data
  #---------------
  Exp <- obs[,2:8]
  y <- obs$Y
  Cov <- obs[,9:12]
  
  
  ### Fit BKMR
  #-----------
  fit_bkmr <- kmbayes(y = y, Z = Exp, X = Cov, iter = 10000, verbose = TRUE, varsel = TRUE)
  

  ### Increase all exposures by 1 unit (g-computation principle)
  #-------------------------------------------------------------
  Z_increased <- Exp + 1
  
  
  ### Compute posterior mean & variance at original and increased exposures
  #------------------------------------------------------------------------
  post_orig <- SamplePred(fit_bkmr, Znew = Exp, method = "exact")
  post_increased <- SamplePred(fit_bkmr, Znew = Z_increased, method = "exact")
  
  
  ### Difference
  #--------------------------------------------------------------
  causal_effect_draws <- rowMeans(post_increased - post_orig)
  
  
  ### Summarize causal effect posterior
  #------------------------------------
  mean_effect <- mean(causal_effect_draws)
  ci_lower <- quantile(causal_effect_draws, 0.025)
  ci_upper <- quantile(causal_effect_draws, 0.975)
  
  
  ### Joint effect summary
  #-----------------------
  ci_width <- ci_upper - ci_lower
  empirical_power <- (ci_lower > 0 | ci_upper < 0)
  CI_coverage <- (ci_lower <= psi[1]) & (ci_upper >= psi[1])
  
  
  ### Store all results
  #--------------------
  parameters <- c("Joint_Increase_All_Exposures", ExtractPIPs(fit_bkmr)[,1])
  estimates <- c(mean_effect, rep(NA,7))
  ci_widths_all <- c(ci_width, rep(NA,7))
  empirical_powers_all <- c(empirical_power, rep(NA,7))
  CI_coverages_all <- c(CI_coverage, rep(NA,7))
  
  
  ### Store results together
  #-------------------------
  coef_bkmr <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = paste0("Bayesian kernel machine regression"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = c(NA, ExtractPIPs(fit_bkmr)[,2])
  )
  
  
  return(coef_bkmr)
}
