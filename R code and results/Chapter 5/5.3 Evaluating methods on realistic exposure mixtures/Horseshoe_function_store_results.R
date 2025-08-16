################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from Horseshoe regression
# 19/06/2025
################################################################################


### Load libraries 
#-----------------
library(nimble)
library(coda)


### Function to store output
#---------------------------
run_Horseshoe <- function(obs, alpha = 0.05, sim_id = 1, psi){
  
tryCatch({ 
  
  ### Define Likelihood & priors
  #-----------------------------
  horseshoeCode  <- nimbleCode({
    
    # Likelihood 
    for (i in 1:n) {
      mu[i] <- beta0 + x[i, 1] * beta[1] + x[i, 2] * beta[2] + x[i, 3] * beta[3] + x[i, 4] * beta[4] + 
        x[i, 5] * beta[5] + x[i, 6] * beta[6] + x[i, 7] * beta[7] + x[i, 8] * beta[8] + x[i, 9] * beta[9] + x[i, 10] * beta[10] + x[i, 11] * beta[11] 
      y[i] ~ dnorm(mu[i], sd = sigma)
    }
    
    # Spike & slab prior on regression coefficients
    for (j in 1:p-5) {
      lambda[j] ~ T(dt(0, 1, df = 1) , 0 , Inf) # heavy-tailed half-Cauchy
      beta[j] ~ dnorm(0, sd = tau * lambda[j])
    }
    
    # Prior on sigma (improper Jeffreys prior)
    log_sigma ~ dnorm(0, sd = 100)  # Approximates π(σ) ∝ 1/σ via transformation
    sigma <- exp(log_sigma)  # Ensure sigma > 0
    
    # Priors
    beta0 ~ dnorm(0, sd = 100)
    beta[8] ~ dnorm(0, sd = 100)
    beta[9] ~ dnorm(0, sd = 100)
    beta[10] ~ dnorm(0, sd = 100)
    beta[11] ~ dnorm(0, sd = 100)
    
    # Hyperpriors
    tau ~ T(dt(0, 1, df = 1) , 0, Inf)
    
    # Joint effect
    sum_beta <- sum(beta[1:7])
  })
  
  
  ### Define constants
  #-------------------
  n <- nrow(obs)
  p <- ncol(obs)
  my.constants <- list(n = n, p = p)
  
  
  ### Initial values
  #-----------------
  my.inits <-  list(
    list(beta0 = 0, beta = rep(0, p-1), log_sigma = log(1), lambda = rep(0.5, 7), tau = 1 ),
    list(beta0 = 0.5, beta = rep(-1, p-1), log_sigma = log(2), lambda = rep(0.5, 7), tau = 1)
  )
  
  
  ### Specify parameters to monitor
  #--------------------------------
  parameters <- c("beta0","beta","sigma","lambda", "sum_beta","tau")
  
  
  ### Run the MCMC chain
  #---------------------
  my.data <- list(
    y = obs$Y,
    x = as.matrix(obs[-1])
  )
  model.sim_horseshoe <- nimbleMCMC(code = horseshoeCode ,
                                    data = my.data,
                                    constants = my.constants,
                                    inits = my.inits,
                                    monitors = parameters,
                                    niter = 100000,
                                    nburnin = 5000,
                                    nchains = 2,
                                    thin = 10,
                                    summary = TRUE,
                                    samplesAsCodaMCMC = TRUE,
                                    WAIC = TRUE)
  
  
  ### Store all coefficients
  #-------------------------
  coefs <- as.data.frame(round(model.sim_horseshoe$summary$all.chains,3))
  rownames(coefs) <-c(names(obs)[-1], "Intercept", "lambda1", "lambda2", "lambda3",
                      "lambda4", "lambda5", "lambda6", "lambda7", 
                      "sigma","sum_beta","tau")
  
  
  ### Joint effect summary
  #-----------------------
  joint_effect <- coefs[21,1]
  empirical_power <- as.numeric(coefs[21, 4] > 0 | coefs[21, 5] < 0)
  CI_coverage <- as.numeric(coefs[21, 4] <= psi[1] & coefs[21, 5] >= psi[1])
  ci_width <- as.numeric(coefs[21, 5] - coefs[21, 4])
  
  
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
    coef_w <- coefs[parameters[j + 1],]
    chemical_estimates[j] <- as.numeric(coef_w[1])
    
    ci <- as.numeric(coef_w[4:5])
    ci_lowers[j] <- ci[1]
    ci_uppers[j] <- ci[2]
    
    ci_widths[j] <- ci[2] - ci[1]
    empirical_powers[j] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[j] <- (ci[1] <= psi[1] * psi[j+1]) & (ci[2] >= psi[1] * psi[j+1])
  }
  
  
  ### Add joint effect to the chemical list
  #----------------------------------------
  estimates <- c(joint_effect, chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Store and return as a dataframe
  #----------------------------------
  coef_horseshoe <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = "Bayesian horseshoe regression",
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = rep(NA,8)
  )
  
  
  return(coef_horseshoe)}, 
  error = function(e) {
    message("Error in run_Horseshoe (sim_id ", sim_id, "): ", e$message)
    return(NULL)
  })
}