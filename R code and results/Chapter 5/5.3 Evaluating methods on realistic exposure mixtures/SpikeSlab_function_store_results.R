################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Function to store all the results from spike and slab
# 05/07/2025
################################################################################


### Load libraries 
#-----------------
library(nimble)


### Function to store output
#---------------------------
run_SpikeSlab <- function(obs, alpha = 0.05, sim_id = 1, psi){
  
  
  ### Define Likelihood & priors
  #-----------------------------
  spikeslabCode <- nimbleCode({
    
    # Likelihood 
    for (i in 1:n) {
      mu[i] <- beta0 + x[i, 1] * beta[1] + x[i, 2] * beta[2] + x[i, 3] * beta[3] + x[i, 4] * beta[4] + 
        x[i, 5] * beta[5] + x[i, 6] * beta[6] + x[i, 7] * beta[7] + x[i, 8] * beta[8] + x[i, 9] * beta[9] + x[i, 10] * beta[10] + x[i, 11] * beta[11] 
      y[i] ~ dnorm(mu[i], sd = sigma)
    }
    
    # Spike & slab prior on regression coefficients
    for (j in 1:(p-5)) {
      beta[j] ~  dnorm(0, sd = lambda[j] * c + (1-lambda[j]) * tau)
      lambda[j] ~ dbern(pi)
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
    
    # Hyperparameters
    pi ~ dunif(0,1)
    tau ~ dgamma(5, 1000)  # mean = 0.005
    c ~ dgamma(2, 20)      # mean = 0.1
    
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
    list(beta0 = 0, beta = rep(0, p-1), log_sigma = log(1), lambda = rep(0, 7), pi =0.8, c = 0.1, tau = 100 ),
    list(beta0 = 0.5, beta = rep(-1, p-1), log_sigma = log(2), lambda = rep(0, 7), pi =0.9, c = 0.1, tau = 100)
  )
  
  
  ### Specify parameters to monitor
  #--------------------------------
  parameters <- c("beta0","beta","sigma","lambda", "sum_beta","pi", "c", "tau")
  
  
  ### Run the MCMC chain
  #---------------------
  my.data <- list(
    y = obs$Y,
    x = as.matrix(obs[-1])
  )
  model.sim_spikeslab <- nimbleMCMC(code = spikeslabCode,
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
  coefs <- as.data.frame(round(model.sim_spikeslab$summary$all.chains,3))
  rownames(coefs) <-c(names(obs)[-1], "Intercept","c", paste0("lambda", names(obs)[2:8]), "pi", 
                      "sigma","sum_beta","tau")
  
  
  ### Joint effect summary
  #-----------------------
  joint_effect <- coefs["sum_beta", 1]
  empirical_power <- as.numeric(coefs["sum_beta", 4] > 0 | coefs["sum_beta", 5] < 0)
  CI_coverage <- as.numeric(coefs["sum_beta", 4] <= psi[1] & coefs["sum_beta", 5] >= psi[1])
  ci_width <- as.numeric(coefs["sum_beta", 5] - coefs["sum_beta", 4])
  
  
  ### Placeholder vectors
  #----------------------
  chemical_estimates <- numeric(7)
  ci_lowers <- numeric(7)
  ci_uppers <- numeric(7)
  ci_widths <- numeric(7)
  pips <- numeric(7)
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
    pips[j] <- as.numeric(coefs[paste0("lambda",parameters[j + 1]),][1])
  }
  
  
  ### Add joint effect to the chemical list
  #----------------------------------------
  estimates <- c(joint_effect, chemical_estimates)
  ci_widths_all <- c(ci_width, ci_widths)
  empirical_powers_all <- c(empirical_power, empirical_powers)
  CI_coverages_all <- c(CI_coverage, CI_coverages)
  
  
  ### Compute Posterior Inclusion Probabilities
  #--------------------------------------------
  pips_all <- c(NA, pips)  # NA for joint effect
  
  
  coef_spikeslab <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = "Bayesian spike & slab regression",
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    PIP = pips_all
  )
  
  
  return(coef_spikeslab)
}

