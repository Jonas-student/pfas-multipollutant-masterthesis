################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from random subset WQS regression
# 24/06/2025
################################################################################


### Load libraries 
#-----------------
library(gWQS)


### Function to store output
#---------------------------
run_WQS_rs <- function(obs, bs = 100, signal = "one" , alpha = 0.05, sim_id = 1, psi, q = NULL){
  
  library(stats)
  
  ### Fit a WQS regression model
  #-----------------------------
  resultwqs <- tryCatch({
    suppressWarnings(
      resultwqs <- gwqs(
        Y ~ wqs + age + gender + Educ2 + Educ3,
        mix_name = c("lbpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
        data = obs,
        q = q,
        rh = 1,
        b = bs,
        rs = TRUE,
        b1_pos = FALSE,
        b_constr = FALSE,
        family = gaussian,
        signal = signal,
      )
    )
    resultwqs
    
  }, error = function(e) {
    message("WQS failed on iteration ", sim_id, ": ", e$message)
    NULL
  })
  
  
  ### Immediately return if the model failed
  #-----------------------------------------
  if (is.null(resultwqs)) {
    return(NULL)
  }
  
  
  ### Joint effect summary
  #-----------------------
  coefs <- summary(resultwqs)$coefficients
  joint_effect <- coefs[2,1]
  empirical_power <- as.numeric(confint(resultwqs)[2, "2.5 %"] > 0 | confint(resultwqs)[2, "97.5 %"] < 0)
  CI_coverage <- as.numeric(confint(resultwqs)[2, "2.5 %"] <= psi[1] & confint(resultwqs)[2, "97.5 %"] >= psi[1])
  ci_width <- as.numeric(confint(resultwqs)[2, "97.5 %"] - confint(resultwqs)[2, "2.5 %"])
  
  
  ### Placeholder vectors
  #----------------------
  chemical_estimates <- numeric(7)
  ci_lowers <- numeric(7)
  ci_uppers <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  parameters <- c("Joint_Increase_All_Exposures", "lbpfos", 
                  "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
  
  
  ### Calculate statistics for each chemical
  #-----------------------------------------
  for(j in seq(1,7,1)){
    chemical_estimates[j] <- resultwqs$final_weights[parameters[j + 1],2]
    
    
    # Cannot be calculated for single partitioning
    ci_widths[j] <- NA
    empirical_powers[j] <- NA
    CI_coverages[j] <- NA
  }
  
  
  ### Add joint effect to the chemical list
  #----------------------------------------
  estimates <- c(joint_effect, chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Store and return as a dataframe
  #----------------------------------
  coef_WQS <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = paste0("WQS (",signal, ") regression (", bs, " random subsets)"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = rep(NA,8)
  )
  
  
  return(coef_WQS)
}
