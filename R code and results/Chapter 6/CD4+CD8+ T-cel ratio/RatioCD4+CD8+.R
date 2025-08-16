################################################################################
# Case study - Jonas Meijerink
# Log(ratio CD4+/CD8+)
# 25/07/2025
################################################################################


rm(list = ls())


### Load libraries 
#-----------------
library(readxl)
library(dplyr)
library(gWQS)
library(tidyverse)
library(MASS)
library(foreach)
library(randomForest)
library(bkmr)
library(BAS)
library(glmnet)
library(nimble)
library(ggplot2)
library(patchwork)
library(coda)
library(MCMCvis)
library(caret)


### Load data
#------------
df_original <- read_excel("~/Dataset_HBM_3M.xlsx")


################################################################################
#---------------------------- Data manipulations ------------------------------#
################################################################################


### Filter on exclusion parameters
#---------------------------------
df <- df_original %>%
  filter(RB_D1==0 & GENEESMIDDEL_SCHILDKLIERAAND==0 & GENEESMIDDEL_DIABETES==0 & 
           GENEESMIDDEL_NIERZIEKTE==0)
# CORTICOSTEROID had only missing values after excluding previous ones


### Make gender a dummy variable (0 jongen/1 meisje)
#---------------------------------------------------
df$GESL <- (df$GESL-1)


### Create dummy variables for age (1=[12.5, 14.5] jaar;2=]14.5, 15.5] jaar;3=> 15.5 jaar)
#-----------------------------------------------------------------------------------------
df$Age2 <- ifelse(df$AGE==2,1,0)
df$Age3 <- ifelse(df$AGE==3,1,0)


### Create dummy variables for income 1=heel erg moeilijk tot moeilijk;
# 2=lukt om rond te komen;3=lukt om comfortabel te leven
#----------------------------------------------------------------------
df$Income2 <- ifelse(df$RONDKOMEN2==2,1,0)
df$Income3 <- ifelse(df$RONDKOMEN2==3,1,0)


### Create dummy variables for BMI 1=(ernstig) ondergewicht;
# 2=normaal gewicht;3=(ernstig) overgewicht 
#-----------------------------------------------------------
df$BMI2 <- ifelse(df$BMIKLAS==2,1,0)
df$BMI3 <- ifelse(df$BMIKLAS==3,1,0)


### DUmmy geboortegewicht 0=>:2.5 KG;1=<2.5 kg (te licht)
# "birthweight"=df$GEBOORTEGEW_LICHT (weggelaten veel missing data)
### Passief roken thuis  of elders (tabak ) of zelf roken? 0=neen;1=ja


### Create data frame
#--------------------
Data <- data.frame("bpfos" = log(df$bpfos_imp) , "lbpfhxs" = log(df$lbpfhxs_imp), 
                   "lbpfoa"=log(df$lbpfoa_imp), "pfba" = log(df$pfba_imp), 
                   "pfda" = log(df$pfda_imp) , "pfna" = log(df$pfna_imp) , 
                   "pfos" = log(df$pfos_imp), "age2"=df$Age2, "age3"=df$Age3,
                   "gender"=df$GESL, "income2"=df$Income2,
                   "income3"=df$Income3 , "bmi2"=df$BMI2, "bmi3"=df$BMI3,
                   "smoke"=df$ALLSMOKE, "birthweight"=df$GEBOORTEGEW_LICHT)
exposures <- data.frame("bpfos" = (df$bpfos_imp) , "lbpfhxs" = (df$lbpfhxs_imp), 
                        "lbpfoa"=(df$lbpfoa_imp), "pfba" = (df$pfba_imp), 
                        "pfda" = (df$pfda_imp) , "pfna" = (df$pfna_imp) , 
                        "pfos" = (df$pfos_imp))
exposures <- na.omit(exposures)


### Extract SDs used for scaling
#-------------------------------
scaled_data <- scale(Data[1:7])
sds <- attr(scaled_data, "scaled:scale")
mu <- attr(scaled_data, "scaled:center")


### Extract Q1 and Q3 for all exposures
#--------------------------------------
Q1 <- apply(exposures[, 1:7], 2, quantile, probs = 0.25)
Q3 <- apply(exposures[, 1:7], 2, quantile, probs = 0.75)


### Center and scale exposures
#-----------------------------
Data[1:7] <- scale(Data[1:7])


### Define the outcome
#---------------------
Y <- log(df$cd4plus_imp/df$cd8plus_imp)


### Delete missing values
#------------------------
obs <- data.frame(Y , Data)
obs <- na.omit(obs)


################################################################################
##
#----------------------- Individual effect estimation -------------------------#
##
################################################################################


### Define outcome, exposures and covariates
#-------------------------------------------
outcome_name <- names(obs)[1]
exposures <- names(obs)[2:8]
covariates <- names(obs)[9:length(obs)]


### Define additive linear model formula
#---------------------------------------
formula_rf <- as.formula(paste(outcome_name, "~", paste(c(exposures, covariates), collapse = " + ")))


#------------------------------------------------------------------------------#
########################### Ridge regression ################################### 
#------------------------------------------------------------------------------#


### Define number of bootstrap samples
#-------------------------------------
n_boot <- 2000


### Prepare data
#---------------
X <- model.matrix(formula_rf, data = obs)[, -1]
Y <- obs$Y
penalty_factors <- c(rep(1, 7), rep(0, 9))


### Plot the coefficients path
#-----------------------------
lambda <- exp(seq(log(1e-6), log(1e3), length.out = 10000))
ridge_fit <- glmnet(X, Y, alpha = 0, penalty.factor = penalty_factors,
                    type.measure = "mse", nfolds = 5, lambda = lambda)
coef_mat <- as.matrix(ridge_fit$beta)[1:7, ]  # exclude intercept
coef_mat <- coef_mat[, ncol(coef_mat):1]
matplot(log10(lambda), t(coef_mat), type = "l", lty = 1, col = 1:7, 
        xlab = "log10(Lambda)", ylab = "Coefficient", main = "Ridge coefficients")
abline(v = -2, lty = 2)
abline(v = -4, lty = 2)
legend("bottomright", legend = c("PFOS (branched)", "PFHXS (total)", "PFOA (total)", "PFBA",
                              "PFDA", "PFNA", "PFOS"), col = 1:7, lty = 1, cex = 1.1, bty = "n")


### Fit Ridge regression using cross-validated lambda
#----------------------------------------------------
cv_fit <- cv.glmnet(X, Y, alpha = 0, penalty.factor = penalty_factors,
                    type.measure = "mse", nfolds = 5, lambda = 
                      exp(seq(log(1e-4), log(1e-1), length.out = 1000)))
ridge_fit <- glmnet(X, Y, alpha = 0, penalty.factor = penalty_factors,
                    lambda = cv_fit$lambda.1se)


### Check model fit and assumptions
#----------------------------------
y_pred <- predict(ridge_fit, newx = rbind(X))
residuals <- Y - y_pred
r_squared <- 1 - sum((Y - y_pred)^2) / sum((Y - mean(Y))^2)
cat("R-squared:", r_squared, "\n")


### Q-Q plot
#-----------
qqnorm(residuals); qqline(residuals, col = "red")


for(i in 2:8){
  par(mfrow=c(2,1))
  ### Homoskedasticity: Residuals vs. Fitted
  #-----------------------------------------
  plot(obs[,i], I(residuals^2), 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Homoskedasticity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], residuals), col = "blue",
        lwd = 2, lty = 1)
  
  
  ### Linearity: Residuals vs. Predictors
  #--------------------------------------
  plot(obs[,i], residuals, 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Linearity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], residuals), col = "blue",
        lwd = 2, lty = 1)
}


### Store bootstrap results
#--------------------------
chemical_effects <- matrix(nrow = n_boot, ncol = 7)


### Run bootstrap sampling
#-------------------------
set.seed(485721023)
for(i in 1:n_boot){
  
  
  # Bootstrap sample
  boot_indices <- sample(1:nrow(obs), replace = TRUE)
  X_boot <- X[boot_indices, ]
  Y_boot <- Y[boot_indices]
  
  
  # Fit Ridge regression using cross-validated lambda
  cv_fit <- cv.glmnet(X_boot, Y_boot, alpha = 0, penalty.factor = penalty_factors,
                      type.measure = "mse", nfolds = 5, lambda = 
                        exp(seq(log(1e-4), log(1e-2), length.out = 1000)))
  ridge_fit <- glmnet(X_boot, Y_boot, alpha = 0, penalty.factor = penalty_factors,
                      lambda = cv_fit$lambda.1se)
  
  
  # Extract coefficients
  coef_ridge <- as.vector(coef(ridge_fit))
  
  
  # Extract exposure estimates and unscale
  beta_unscaled <- coef_ridge[2:8] /sds
  chemical_effects[i,] <- (Q3/Q1)**beta_unscaled
}


### Parameter names
#------------------
chemicals <- c("bpfos", "lbpfhxs ", "lbpfoa ", "pfba ", "pfda ", "pfna ", "pfos")


### Placeholders for chemical results
#------------------------------------
chemical_estimates <- numeric(7)
ci <- cbind(numeric(7),numeric(7))


### Calculate chemical-specific statistics
#-----------------------------------------
for (i in 1:7) {
  effect_samples <- chemical_effects[, i]
  chemical_estimates[i] <- mean(effect_samples)
  ci[i,] <- quantile(effect_samples, probs = c(0.025, 0.975))
}
coefs_ridge <- data.frame(
  "parameter" = chemicals,
  "Mean" = chemical_estimates,
  "95%CI_low" = ci[,1],
  "95%CI_upp" = ci[,2],
  "method"="Ridge regression"
)


#------------------------------------------------------------------------------#
############################# OLS regression ###################################
#------------------------------------------------------------------------------#
# Model assumptions are checked bellow for OLS for joint effect estimation


### Store bootstrap results
#--------------------------
chemical_effects <- matrix(nrow = n_boot, ncol = 7)


### Run bootstrap sampling
#-------------------------
set.seed(485721023)
for(i in 1:n_boot){
  
  
  # Bootstrap sample
  boot_indices <- sample(1:nrow(obs), replace = TRUE)
  data_boot <- obs[boot_indices,]
  
  
  # Fit OLS regression
  model <- lm(formula_rf, data = data_boot)
  
  
  # Extract coefficients
  coef_ols <- summary(model)$coefficients[2:8]
  
  
  # Extract exposure estimates and unscale
  beta_unscaled <- coef_ols /sds
  chemical_effects[i,] <- (Q3/Q1)**beta_unscaled
}


### Parameter names
#------------------
chemicals <- c("bpfos", "lbpfhxs ", "lbpfoa ", "pfba ", "pfda ", "pfna ", "pfos")


### Placeholders for chemical results
#------------------------------------
chemical_estimates <- numeric(7)
ci <- cbind(numeric(7),numeric(7))


### Calculate chemical-specific statistics
#-----------------------------------------
for (i in 1:7) {
  effect_samples <- chemical_effects[, i]
  chemical_estimates[i] <- mean(effect_samples)
  ci[i,] <- quantile(effect_samples, probs = c(0.025, 0.975))
}
coefs_ols <- data.frame(
  "parameter" = chemicals,
  "Mean" = chemical_estimates,
  "95%CI_low" = ci[,1],
  "95%CI_upp" = ci[,2],
  "method"="Multiple pollutant linear model (OLS)"
)


#------------------------------------------------------------------------------#
###################### Bayesian (flat prior) regression ########################
#------------------------------------------------------------------------------#


### Define Likelihood & priors
#-----------------------------
BayesianCode  <- nimbleCode({
  
  # Likelihood 
  for (i in 1:n) {
    mu[i] <- beta0 + x[i, 1] * beta[1] + x[i, 2] * beta[2] + x[i, 3] * beta[3] + x[i, 4] * beta[4] + 
      x[i, 5] * beta[5] + x[i, 6] * beta[6] + x[i, 7] * beta[7] + x[i, 8] * beta[8] + x[i, 9] * beta[9] + 
      x[i, 10] * beta[10] + x[i, 11] * beta[11] + x[i, 12] * beta[12] + x[i, 13] * beta[13] + 
      x[i, 14] * beta[14] + x[i, 15] * beta[15]+ x[i, 16] * beta[16] 
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
  
  # Spike & slab prior on regression coefficients
  for (j in 1:7) {
    beta[j] ~ dnorm(0, sd = 100)
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
  beta[12] ~ dnorm(0, sd = 100)
  beta[13] ~ dnorm(0, sd = 100)
  beta[14] ~ dnorm(0, sd = 100)
  beta[15] ~ dnorm(0, sd = 100)
  beta[16] ~ dnorm(0, sd = 100)
  
  # Posterior samples
  for(i in 1:n) {
    y_rep[i] ~ dnorm(mu[i], sd = sigma)
  }
})


### Define constants
#-------------------
n <- nrow(obs)
p <- ncol(obs)
my.constants <- list(n = n, p = p)


### Initial values
#-----------------
my.inits <-  list(
  list(beta0 = 0, beta = rep(0, p-1), log_sigma = log(1)),
  list(beta0 = 0.5, beta = rep(-1, p-1), log_sigma = log(2))
)


### Specify parameters to monitor
#--------------------------------
parameters <- c("beta0","beta","sigma","mu","y_rep")


### Run the MCMC chain
#---------------------
my.data <- list(
  y = obs$Y,
  x = as.matrix(obs[-1])
)
set.seed(485721023)
model.sim_bayesian <- nimbleMCMC(code = BayesianCode ,
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


### Check convergence
#--------------------
mcmc_list <- coda::mcmc.list(
  lapply(model.sim_bayesian$samples, mcmc)
)
plot(mcmc_list[,c(1:17)]) 
gelman.diag(mcmc_list[,c(1:17)])
effectiveSize(mcmc_list[,c(1:17)])


### Assumption Checks
#-------------------------------------------------------------------------------
mu <- MCMCpstr(mcmc_list, params = "mu")$mu
resid <- obs$Y - mu


### Normality assumption
#-----------------------
qqnorm(resid, main = "Q-Q Plot for Residuals")
qqline(resid, col = "red")


for(i in 2:8){
  par(mfrow=c(2,1))
  ### Homoskedasticity: Residuals vs. Fitted
  #-----------------------------------------
  plot(obs[,i], I(resid^2), 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Homoskedasticity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], resid), col = "blue",
        lwd = 2, lty = 1)
  
  
  ### Linearity: Residuals vs. Predictors
  #--------------------------------------
  plot(obs[,i], resid, 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Linearity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], resid), col = "blue",
        lwd = 2, lty = 1)
}


### Density Overlay
#------------------
y_rep <- rbind(mcmc_list$chain1[,278:536],mcmc_list$chain2[,278:536])
y_rep_long <- data.frame(Value = as.vector(t(y_rep)), Type = "Posterior predictive distribution")
obs_long <- data.frame(Value = obs$Y, Type = "Observed distribution")  # Use log_y
combined_df <- rbind(obs_long, y_rep_long)
ggplot(combined_df, aes(x = Value, color = Type, linetype = Type)) +
  geom_density(size = 1) +
  scale_color_manual(values = c("black", "red")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(title = "Observed vs. posterior predicted distribution of log(y)", x = "log(y)", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


### Posterior predictive p-value
#-------------------------------
obs_mean <- mean(obs$Y)
pred_means <- rowMeans(y_rep)
p_value <- mean(pred_means > obs_mean)
cat("Posterior predictive p-value for mean:", p_value, "\n")


### Compute fold change
#----------------------
samples_matrix <- as.matrix(model.sim_bayesian$samples)
beta_cols <- grep("^beta\\[[1-7]\\]$", colnames(samples_matrix))
beta_samples <- samples_matrix[, beta_cols]
beta_unscaled <- sweep(beta_samples, 2, sds[1:7], FUN = "/")
fold_changes <- (Q3/Q1)**beta_unscaled


### Store all coefficients
#-------------------------
coefs_bayes <- data.frame(
  "parameter" = chemicals,
  "Mean" = apply(fold_changes, 2, mean),
  "95%CI_low" = apply(fold_changes, 2, quantile, probs = 0.025),
  "95%CI_upp" = apply(fold_changes, 2, quantile, probs = 0.975),
  "method" = "Bayesian regression (flat prior)"
)


#------------------------------------------------------------------------------#
############################ Horseshoe regression ##############################
#------------------------------------------------------------------------------#


### Define Likelihood & priors
#-----------------------------
horseshoeCode  <- nimbleCode({
  
  # Likelihood 
  for (i in 1:n) {
    mu[i] <- beta0 + x[i, 1] * beta[1] + x[i, 2] * beta[2] + x[i, 3] * beta[3] + x[i, 4] * beta[4] + 
      x[i, 5] * beta[5] + x[i, 6] * beta[6] + x[i, 7] * beta[7] + x[i, 8] * beta[8] + x[i, 9] * beta[9] + 
      x[i, 10] * beta[10] + x[i, 11] * beta[11] + x[i, 12] * beta[12] + x[i, 13] * beta[13] + 
      x[i, 14] * beta[14] + x[i, 15] * beta[15]+ x[i, 16] * beta[16]
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
  
  # Spike & slab prior on regression coefficients
  for (j in 1:7) {
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
  beta[12] ~ dnorm(0, sd = 100)
  beta[13] ~ dnorm(0, sd = 100)
  beta[14] ~ dnorm(0, sd = 100)
  beta[15] ~ dnorm(0, sd = 100)
  beta[16] ~ dnorm(0, sd = 100)
  
  # Hyperpriors
  tau ~ T(dt(0, 1, df = 1) , 0, Inf)
  
  
  # Posterior predictive distribution
  for(i in 1:n) {
    y_rep[i] ~ dnorm(mu[i], sd = sigma)
  }
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
parameters <- c("beta0","beta","sigma","lambda","tau","mu","y_rep")


### Run the MCMC chain
#---------------------
my.data <- list(
  y = obs$Y,
  x = as.matrix(obs[-1])
)
set.seed(485721023)
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


### Check convergence
#--------------------
mcmc_list <- mcmc.list(
  lapply(model.sim_horseshoe$samples, mcmc)
)
plot(mcmc_list[,c(1:24,285:286)]) 
gelman.diag(mcmc_list[,c(1:24,285:286)])
effectiveSize(mcmc_list[,c(1:24,285:286)])


### Assumption Checks
#-------------------------------------------------------------------------------
mu <- MCMCpstr(mcmc_list, params = "mu")$mu
resid <- obs$Y - mu


### Normality assumption
#-----------------------
qqnorm(resid, main = "Q-Q Plot for Residuals")
qqline(resid, col = "red")


for(i in 2:8){
  par(mfrow=c(2,1))
  ### Homoskedasticity: Residuals vs. Fitted
  #-----------------------------------------
  plot(obs[,i], I(resid^2), 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Homoskedasticity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], resid), col = "blue",
        lwd = 2, lty = 1)
  
  
  ### Linearity: Residuals vs. Predictors
  #--------------------------------------
  plot(obs[,i], resid, 
       xlab = names(obs)[i], ylab = "Residuals", 
       main = "Linearity Check: log transformation")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(loess.smooth(obs[,i], resid), col = "blue",
        lwd = 2, lty = 1)
}


### Density Overlay
#------------------
y_rep <- cbind(mcmc_list$chain1[,286:544],mcmc_list$chain2[,286:544])
y_rep_long <- data.frame(Value = as.vector(t(y_rep)), Type = "Posterior predictive distribution")
obs_long <- data.frame(Value = obs$Y, Type = "Observed distribution")  # Use log_y
combined_df <- rbind(obs_long, y_rep_long)
ggplot(combined_df, aes(x = Value, color = Type, linetype = Type)) +
  geom_density(size = 1) +
  scale_color_manual(values = c("black", "red")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(title = "Observed vs. posterior predicted distribution of log(y)", x = "log(y)", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


### Posterior predictive p-value
#-------------------------------
obs_mean <- mean(obs$Y)
pred_means <- colMeans(y_rep)
p_value <- mean(pred_means > obs_mean)
cat("Posterior predictive p-value for mean:", p_value, "\n")


### Calcluate R squared
#----------------------
SS_res <- sum((obs$Y - mu)^2)
SS_tot <- sum((obs$Y - mean(obs$Y))^2)
R2 <- 1 - SS_res / SS_tot
print(R2)


### Calculate fold change
#------------------------
samples_matrix <- as.matrix(model.sim_horseshoe$samples)
beta_cols <- grep("^beta\\[[1-7]\\]$", colnames(samples_matrix))
beta_samples <- samples_matrix[, beta_cols]
beta_unscaled <- sweep(beta_samples, 2, sds[1:7], FUN = "/")  # sds should be in same order
fold_changes <- (Q3/Q1)**beta_unscaled


### Store all coefficients
#-------------------------
coefs_horse <- data.frame(
  "parameter" = chemicals,
  "Mean" = apply(fold_changes, 2, mean),
  "95%CI_low" = apply(fold_changes, 2, quantile, probs = 0.025),
  "95%CI_upp" = apply(fold_changes, 2, quantile, probs = 0.975),
  "method"="Bayesian horseshoe regression"
)


#------------------------------------------------------------------------------#
######################### Repeated holdout WQS #################################
#------------------------------------------------------------------------------#


### Fit repeated holdout WQS regression model (negative direction)
#-----------------------------------------------------------------
resultwqs_neg <- gwqs(
  Y ~ wqs + age2 + age3 + gender + income2 + income3 + bmi2 + bmi3 + 
    smoke + birthweight,
  mix_name = c("bpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
  data = obs,
  q = NULL,
  rh = 100,
  b = 100,
  b1_pos = FALSE,
  b_constr = FALSE,
  family = gaussian,
  seed = 485721023
)


### Scatter plot y vs wqs
#------------------------
gwqs_scatterplot(resultwqs_neg)


### Scatter plot residuals vs fitted values
#------------------------------------------
gwqs_fitted_vs_resid(resultwqs_neg)



### Check significance of joint effect
#-------------------------------------
resultwqs_neg$fit$coefficients

# Was not found to be significant so stop the analysis, weights cannot be 
# interpreted.


### Fit repeated holdout WQS regression model (positive direction)
#-----------------------------------------------------------------
resultwqs_pos <- gwqs(
  Y ~ wqs + age2 + age3 + gender + income2 + income3 + bmi2 + bmi3 + 
    smoke + birthweight,
  mix_name = c("bpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
  data = obs,
  q = NULL,
  rh = 100,
  b = 100,
  b1_pos = TRUE,
  b_constr = TRUE,
  family = gaussian,
  seed = 485721023
)


### Scatter plot y vs wqs
#------------------------
gwqs_scatterplot(resultwqs_pos)


### Scatter plot residuals vs fitted values
#------------------------------------------
gwqs_fitted_vs_resid(resultwqs_pos)



### Check significance of joint effect
#-------------------------------------
resultwqs_pos$fit$coefficients

# Was not found to be significant so stop the analysis, weights cannot be 
# interpreted.


#------------------------------------------------------------------------------#
################### Bayesian kernel machine regression #########################
#------------------------------------------------------------------------------#  


### Prepare data
#---------------
Exp <- obs[,2:8]
y <- obs$Y
Cov <- obs[,9:17]


### Fit BKMR
#-----------
set.seed(485721023)
fit_bkmr <- kmbayes(y = y, Z = Exp, X = Cov, iter = 10000, verbose = TRUE, 
                    varsel = TRUE, control.params = list(
                      a.p0 = 1,
                      b.p0 = 1
                    ))


### Investigate convergence
#--------------------------
for(i in 1:9){
  TracePlot(fit = fit_bkmr, par = "beta", comp = i)
}
for(i in 1:7){
  TracePlot(fit = fit_bkmr, par = "r", comp = i)
}
TracePlot(fit = fit_bkmr, par = "sigsq.eps")


### Extract posterior inclusion probabilities
#--------------------------------------------
ExtractPIPs(fit_bkmr)


################################################################################
#-------------------------- Joint effect estimation ---------------------------#
################################################################################


#------------------------------------------------------------------------------#
############################### Random forest ##################################
#------------------------------------------------------------------------------#


### Define parameters
#--------------------
mtry_values <- 2:6
ntree_vals <- seq(10, 500, by = 10)
k_folds <- 5
n_repeats <- 20
results <- list()
count <- 1


### Explore hyperparamer choices in terms of MSE
#-----------------------------------------------
set.seed(1234)
for (repeat_i in 1:n_repeats) {
  
  
  # Create folds for this repeat
  folds <- createFolds(obs$Y, k = k_folds, list = TRUE, returnTrain = FALSE)
  
  for (mtry_val in mtry_values) {
    for (ntree_val in ntree_vals) {
      fold_mse <- numeric(k_folds)
      
      for (fold_i in seq_along(folds)) {
        train_idx <- setdiff(1:nrow(obs), folds[[fold_i]])
        test_idx <- folds[[fold_i]]
        
        train_data <- obs[train_idx, ]
        test_data <- obs[test_idx, ]
        
        rf <- randomForest(formula_rf, data = train_data,
                           mtry = mtry_val, ntree = ntree_val)
        
        pred <- predict(rf, newdata = test_data)
        fold_mse[fold_i] <- mean((pred - test_data$Y)^2)
      }
      
      results[[count]] <- data.frame(
        "repeat" = repeat_i,
        mtry = mtry_val,
        ntree = ntree_val,
        mean_fold_mse = mean(fold_mse)
      )
      count <- count + 1
    }
  }
}

### Summarize over repeats
#-------------------------
df_results <- bind_rows(results)
df_summary <- df_results %>%
  group_by(mtry, ntree) %>%
  summarise(
    mean_cv_mse = mean(mean_fold_mse),
    se_cv_mse = sd(mean_fold_mse) / sqrt(n_repeats),
    .groups = "drop"
  )


### Plot
#-------
ggplot(df_summary, aes(x = ntree, y = mean_cv_mse, color = factor(mtry))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_cv_mse - se_cv_mse,
                  ymax = mean_cv_mse + se_cv_mse,
                  fill = factor(mtry)), alpha = 0.2, color = NA) +
  labs(
    title = paste0(n_repeats, "-Repeated ", k_folds, "-Fold CV Test MSE"),
    subtitle = "Random Forest, varying mtry",
    x = "Number of Trees",
    y = "Mean CV MSE",
    color = "mtry",
    fill = "mtry"
  ) +
  theme_minimal(base_size = 14)


### Step 1: Create counterfactual datasets
#-----------------------------------------
cf_low  <- obs[-1]
cf_high <- obs[-1]


### Step 2: Replace only the exposures (cols 1:7) with Q1 or Q3
#--------------------------------------------------------------
mu <- attr(scaled_data, "scaled:center")
for (j in 1:259) {
  cf_low[j,1:7]  <- (log(Q1)-mu)/sds
  cf_high[j,1:7] <- (log(Q3)-mu)/sds
}


### Function to compute ATE for a given data set
#-----------------------------------------------
compute_ate <- function(data) {
  
  
  ### Step 3: Fit Random Forest model with exposures + covariates
  #--------------------------------------------------------------
  rf_model <- randomForest(formula_rf, data = data, ntree = 500, mtry=2)
  
  
  ### Step 5: Predict using full model (with covariates unchanged)
  #---------------------------------------------------------------
  pred_low  <- predict(rf_model, newdata = cf_low)
  pred_high <- predict(rf_model, newdata = cf_high)
  
  
  ### Step 6: Compute average causal effect
  #----------------------------------------
  joint_effect <- mean(exp(pred_high))/mean(exp(pred_low))
  
  
  return(joint_effect)
}


### Vector for bootstrap results
#-------------------------------
n_boot <- 2000
ate_boot <- numeric(n_boot)


### Bootstrapping
#----------------
set.seed(485721023)
for (b in 1:n_boot) {
  boot_indices <- sample(1:nrow(obs), nrow(obs), replace = TRUE)
  boot_sample <- obs[boot_indices, ]
  ate_boot[b] <- compute_ate(boot_sample)
}


### Calculate bootstrap statistics
#---------------------------------
alpha<-0.05
ci_lower <- quantile(ate_boot, alpha / 2)
ci_upper <- quantile(ate_boot, 1 - alpha / 2)
ate <- mean(ate_boot)
coefs_Gcomp <- data.frame(
  "parameter" = "Joint effect",
  "Mean" = ate,
  "95%CI_low" = ci_lower,
  "95%CI_upp" = ci_upper,
  "method" = "Random forest g-computation"
)


#------------------------------------------------------------------------------#
########################## OLS sum of joint effects ############################
#------------------------------------------------------------------------------#


###-### Step 1: fit a simple linear additive regression model
#-------------------------------------------------------------------------------


### Fit linear model
#-------------------
model <- lm(formula_rf, data = obs)


# Extract residuals and fitted values
#-------------------------------------
residuals <- resid(model)
fitted <- fitted(model)


par(mfrow = c(2,2))
# Plot residuals versus bpfos
#-----------------------------
plot(model$model$bpfos, residuals, 
     xlab = "bpfos", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$bpfos, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the lbpfhxs
#----------------------------------

plot(model$model$lbpfhxs, residuals, 
     xlab = "lbpfhxs", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$lbpfhxs, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the lbpfoa
#---------------------------------
plot(model$model$lbpfoa, residuals, 
     xlab = "lbpfoa", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$lbpfoa, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the pfba
#-------------------------------
plot(model$model$pfba, residuals, 
     xlab = "pfba", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$pfba, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the pfda
#-------------------------------
plot(model$model$pfda, residuals, 
     xlab = "pfda", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$pfda, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the lbpfhxs
#-----------------------------------
plot(model$model$pfna, residuals, 
     xlab = "pfna", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$pfna, residuals), col = "blue",
      lwd = 2, lty = 1)


# Plot residuals versus the lbpfhxs
#-----------------------------------
plot(model$model$pfos, residuals, 
     xlab = "pfos", ylab = "Residuals", 
     main = "Linearity Check: log transformation")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(model$model$pfos, residuals), col = "blue",
      lwd = 2, lty = 1)


# Residuals vs. fitted values
#-----------------------------
plot(fitted, residuals, 
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(loess.smooth(fitted, residuals), col = "blue",
      lwd = 2, lty = 1)


par(mfrow = c(2,1))
# Distribution of residuals
#---------------------------
hist(residuals, main = "Histogram of residuals", freq = F, col = "pink",breaks = 12,)


# Normality assumption
#---------------------------
qq = qqPlot(residuals(model), xlab = "Expected quantiles under normality", 
            ylab = "Empirical quantiles", id = TRUE, envelope = 0.99)


###-### Step 3: estimate joint effect
#-------------------------------------------------------------------------------


### Function to compute ATE for a given data set
#-----------------------------------------------
compute_ate_lm <- function(data) {
  
  
  ### Step 4: Fit linear model (with interactions) 
  #-----------------------------------------------
  lm_model <- lm(formula_rf, data = data)
  
  
  ### Step 5: extract coefficients
  #-------------------------------
  coef <- coef(lm_model)[2:8]
  
  
  ### Calculate joint effect
  #-------------------------
  joint_effect <-1
  for(i in 1:7){
    joint_effect <- (Q3[i]/Q1[i])**(coef[i]/sds[i])
  }
  
  
  return(joint_effect)
}


### Vector for bootstrap results
#-------------------------------
n_boot <- 2000
ate_boot <- numeric(n_boot)


### Bootstrapping
#----------------
set.seed(485721023)
for (b in 1:n_boot) {
  boot_indices <- sample(1:nrow(obs), nrow(obs), replace = TRUE)
  boot_sample <- obs[boot_indices, ]
  ate_boot[b] <- compute_ate_lm(boot_sample)
}


### Calculate bootstrap statistics
#---------------------------------
alpha<-0.05
ci_lower <- quantile(ate_boot, alpha / 2)
ci_upper <- quantile(ate_boot, 1 - alpha / 2)
ate <- mean(ate_boot)
coefs_olsjoint <- data.frame(
  "parameter" = "Joint effect",
  "Mean" = ate,
  "95%CI_low" = ci_lower,
  "95%CI_upp" = ci_upper,
  "method" = "Multiple pollutant linear model (OLS)"
)


################################################################################
#                                                                              #
#---------------------------------- Visuals -----------------------------------#
#                                                                              #
################################################################################


### Compare ridge and horseshoe regression
#-----------------------------------------
df_compare <- bind_rows(coefs_bayes[1:7,], coefs_horse[1:7,],coefs_ridge[1:7,],coefs_ols, coefs_olsjoint,coefs_Gcomp)
df_compare$parameter <- factor(df_compare$parameter, levels = unique(df_compare$parameter))
df_compare$method <- factor(df_compare$method, levels = c("Ridge regression", "Multiple pollutant linear model (OLS)","Bayesian horseshoe regression","Bayesian regression (flat prior)", "Random forest g-computation"))

ggplot(df_compare, aes(x = Mean, y = parameter, color = method)) +
  geom_point(position = position_dodge(width = -0.7), size = 3) +
  geom_errorbarh(aes(xmin = `X95.CI_low`, xmax = `X95.CI_upp`), height = 0.3, position = position_dodge(width = -0.7)) +
  labs(x = "Average interquartile fold change with 95% confidence or credible interval", y = NULL, 
       color = "Method:", title="CD4+/CD8+ T-cell ratio") +
  scale_color_manual(values = c(
    "Bayesian horseshoe regression" = "#FF7F0E",  
    "Ridge regression" = "#1F77B4",    
    "Multiple pollutant linear model (OLS)" = "#1F77B460",  
    "Bayesian regression (flat prior)" = "#FF7F0E65",
    "Random forest g-computation" = "purple"
  )) +
  theme_minimal()+
  scale_y_discrete(labels = c("PFOS (branched)", "PFHXS (total)","PFOA (total)","PFBA","PFDA",  "PFNA","PFOS", "Joint effect")) +
  theme(
    plot.title = element_text(size = 17),
    axis.text.y = element_text(size = 16),        # Remove x-axis labels
    axis.title.y = element_text(size = 17),    # Move legend to bottom
    legend.title = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_blank(),
    legend.position = "bottom"
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  geom_vline(xintercept = 1   , linetype = "dashed", color = "grey60")


