################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Outliers + heavy tailed distributions
# 23/06/2025
################################################################################


rm(list = ls())


### Load packages
#----------------
library(dplyr)
library(ggplot2)
library(MASS)
library(tidyverse)
library(mvtnorm)
library(gWQS)
library(patchwork)
library(grid)
library(hrbrthemes)
library(readxl)


########################### Exposure scenario (c) ##############################
#-------------------------------------------------------------------------------
# Exposure profiles (rows) resampled from the case study dataset to preserve 
# joint distributions and including confounders.
#-------------------------------------------------------------------------------


### Function to store output
#---------------------------
run_WQS_rh_conf <- function(obs, rh = 100, bs = 20, signal = "one" , alpha = 0.05, 
                            sim_id = 1, psi, q = NULL, rho){
  
  
  ### Fit a WQS regression model
  #-----------------------------
  resultwqs <- tryCatch({
    suppressWarnings(
    resultwqs <- gwqs(
      Y ~ wqs + age + gender + Educ2 + Educ3,
      mix_name = c("lbpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
      data = obs,
      q = q,
      rh = rh,
      b = bs,
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
  empirical_power <- as.numeric(coefs[2, "2.5 %"] > 0 | coefs[2, "97.5 %"] < 0)
  CI_coverage <- as.numeric(coefs[2, "2.5 %"] <= psi[1] & coefs[2, "97.5 %"] >= psi[1])
  ci_width <- as.numeric(coefs[2, "97.5 %"] - coefs[2, "2.5 %"])
  
  
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
    coef_w <- resultwqs$final_weights[parameters[j + 1],]
    chemical_estimates[j] <- as.numeric(coef_w[2])
    
    ci <- as.numeric(coef_w[3:4])
    ci_lowers[j] <- ci[1]
    ci_uppers[j] <- ci[2]
    
    ci_widths[j] <- ci[2] - ci[1]
    empirical_powers[j] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[j] <- (ci[1] <= psi[j+1]) & (ci[2] >= psi[j+1])
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
    Method = paste0("WQS (q=",q, ") regression (rh=", rh, " , bs=", bs, ")"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all
  )
  
  
  return(coef_WQS)
}


### Load data for a realistic simulation study
#---------------------------------------------
df <- read_excel("~/Dataset_HBM_3M.xlsx")


### Make gender a dummy variable (0/1)
#-------------------------------------
df$GESL <- (df$GESL-1)


### Create data frame
#--------------------
Data <- data.frame("lbpfos" = log(df$lbpfos_imp) , "lbpfhxs" = log(df$lbpfhxs_imp), 
                   "lbpfoa"=log(df$lbpfoa_imp),"pfba" = log(df$pfba_imp), 
                   "pfda" = log(df$pfda_imp) , "pfna" = log(df$pfna_imp) , 
                   "pfos" = log(df$pfos_imp), "age"=df$OB_Leeftijd, 
                   "highest education"=df$ISCED_HH , "gender"=df$GESL)


### Create dummy variables for education
#---------------------------------------
Data$Educ2 <- ifelse(Data$highest.education==2,1,0)
Data$Educ3 <- ifelse(Data$highest.education==3,1,0)
Data$highest.education <- NULL


### Select continuous variables to simulate
#------------------------------------------
PFAS_clean <- na.omit(Data)


### Convert to matrix
#--------------------
data_matrix <- as.matrix(PFAS_clean)


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#

### Set seed
#-----------
set.seed(123)


### Define parameters
#--------------------
sample_sizes <- c(300)
n_sim <- 50
p <- 7
mu <- c(2.1223 , -1.9563 , -1.960, -0.5490 , -1.3725, 0.13153,1.1773 )
alpha_level <- 0.05
psi <- c(-0.2 , 0.186305814, 0.421321419,  0.051494073, 0.085073273, 0.002444214,
           0.232626663, 0.020734543)


### Store results
#----------------
bootstrap <- data.frame()


for (n in sample_sizes) {
    
    
  ### Repeat this procedure x times
  #--------------------------------
  for (i in 1:n_sim) {
      
    print(i)
    
    
    ### Resample indices (bootstrap resampling)
    #------------------------------------------
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    
    
    ### Extract resampled pseudo-observations
    #----------------------------------------
    simulated_data <- data_matrix[boot_indices, ]
    
    
    ### Convert to data frame and name variables
    #-------------------------------------------
    X <- as.data.frame(simulated_data)
    colnames(X) <- colnames(PFAS_clean)
    X[,1:7] <- scale(X[,1:7])
    
    ### Define the mixture variable names (excluding age)
    #----------------------------------------------------
    mix_vars <- c("lbpfos", "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")

    
    ### Sample the error terms using a variance of 0.33 from classical shrinkage methods
    #-----------------------------------------------------------------------------------
    eps <- rnorm(n , mean = 0, sd = sqrt(0.8466327))
    
    
    ### Use a WQS inspired model to simulate the outcome
    #---------------------------------------------------
    WQS <- 0.186305814*X$lbpfos + 0.051494073*X$pfda + 0.421321419* X$pfba + 
    0.020734543*X$pfos + 0.232626663*X$lbpfoa + 0.002444214*X$pfna + 0.085073273*X$lbpfhxs
    
    
    ### Combine as a linear model functional form
    #--------------------------------------------
    Y <- -1.471401 + psi[1] * WQS - 0.016248 * X$age + 0.168078 * X$gender + 
      0.103059 * X$Educ2 + 0.023381 * X$Educ3 + eps
    obs <- data.frame(Y , X)
      
    
    ### Merge results from all models
    #--------------------------------
    coef_WQS_rh <- run_WQS_rh_conf(obs, rh = 100, bs = 20, signal = "one" , 
                              alpha = 0.05, sim_id = i, psi, q = NULL)
    
    coef_WQS_rh_4 <- run_WQS_rh_conf(obs, rh = 100, bs = 20, signal = "one" , 
                                alpha = 0.05, sim_id = i, psi, q = 4)
    
    coef_WQS_rh_10 <- run_WQS_rh_conf(obs, rh = 100, bs = 20, signal = "one" , 
                                 alpha = 0.05, sim_id = i, psi, q = 10)
    bootstrap <- rbind(bootstrap , coef_WQS_rh,coef_WQS_rh_4,coef_WQS_rh_10 )
  }
}


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Switch model names
#---------------------
bootstrap <- bootstrap %>%
  mutate(
    Quantile = case_when(
      Method == "WQS (q=4) regression (rh=100 , bs=20)" ~ "Quartile",
      Method == "WQS (q=10) regression (rh=100 , bs=20)" ~ "Decile",
      Method == "WQS (q=) regression (rh=100 , bs=20)" ~ "Continuous"
    )
  )


### Select joint effect
#----------------------
coefwqs <- subset(bootstrap, parameter == "Joint_Increase_All_Exposures")
coefwqs$Quantile <- factor(coefwqs$Quantile, levels = c("Quartile", "Decile", "Continuous"))


### Calculate empirical power
#----------------------------
power_summary <- coefwqs %>%
  group_by(Method, Size,Quantile) %>%
  summarise(
    Positives = sum(Power == 1, na.rm = TRUE),  # Count of significant results
    n = n(),                                    # Total number of simulations
    Power = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )



### Make first barplot for scenario (c)
#--------------------------------------
p1<-power_summary %>%
  ggplot(aes(x = Quantile, y = Power, fill = Quantile, color = Quantile)) +  # Add color here
  
  geom_bar(stat = "identity", alpha = 0.5, position = "dodge") +
  
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 17),
    axis.title.y = element_text(size = 18)) +
  ylim(0,1)+
  ylab("Empirical power of the joint effect") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ', 
       color = 'Repeated holdout WQS regression with exposure type: ')


### Select all chemicals
#-----------------------
coefPFAS <- subset(bootstrap, parameter == "lbpfos" | parameter == "pfba" |
                 parameter == "pfda" |parameter == "lbpfhxs" |
                parameter == "pfna" |parameter == "lbpfoa"|parameter == "pfos")


### Allow for calculating the bias
#---------------------------------
coefPFAS <- coefPFAS %>%
  mutate(
    # Match the chemical index in the psi vector
    psi_index = case_when(
      parameter == "lbpfos" ~ 2,
      parameter == "pfba" ~ 3,
      parameter == "pfda" ~ 4,
      parameter == "lbpfhxs" ~ 5,
      parameter == "pfna" ~ 6,
      parameter == "lbpfoa" ~ 7,
      parameter == "pfos" ~ 8,
    ),
    TrueEffect = psi[psi_index]
    ,
    AbsoluteBias = Estimate - TrueEffect,
    RelativeBias = (Estimate - TrueEffect) / abs(TrueEffect)
  )


### Change name labels
#---------------------
name_map <- c("lbpfos" = "PFOS (total) (19%)", "pfda"="PFDA (5%)",
              "pfba" = "PFBA (42%)", "pfos" = "PFOS (2%)",
              "lbpfoa" = "PFOA (total) (23%)", "pfna"="PFNA (~0%)",
               "lbpfhxs" = "PFHXS (total) (9%)" )
coefPFAS$parameter<-name_map[coefPFAS$parameter]
coefPFAS$parameter<-factor(coefPFAS$parameter , name_map)


### Second plot of the individual weights
#--------------------------------------
p2 <- ggplot(coefPFAS, aes(x = Quantile, y = AbsoluteBias, fill = Quantile)) +
  geom_boxplot(alpha=0.5, show.legend = FALSE) +
  geom_jitter(aes(color = Quantile), width = 0.2, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
  theme_ipsum() +
  ylim(-0.4,0.4)+
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 8),
    axis.title.y = element_text(size = 18)
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
  ylab("Absolute bias") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ',
       color = 'Repeated holdout WQS regression with exposure type: ') +
  geom_hline(yintercept = 0, col = "red")


########################### Exposure scenario (a) ##############################
#-------------------------------------------------------------------------------
# Samples drawn from a multivariate t-distribution with two degrees of freedom, char-
# acterised by heavy tails with either no correlation or high correlation among the
# exposures.
#-------------------------------------------------------------------------------


### Function to store output
#---------------------------
run_WQS_rh <- function(obs, rh = 100, bs = 20, signal = "one" , alpha = 0.05, sim_id = 1, psi, q = NULL, rho){
  
  
  ### Fit a WQS regression model
  #-----------------------------
  resultwqs <- tryCatch({
    suppressWarnings(
    resultwqs <- gwqs(
      Y ~ wqs ,
      mix_name = c("X1","X3", "X2", "X7", "X6", "X5", "X4"),
      data = obs,
      q = q,
      rh = rh,
      b = bs,
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
  empirical_power <- as.numeric(coefs[2, "2.5 %"] > 0 | coefs[2, "97.5 %"] < 0)
  CI_coverage <- as.numeric(coefs[2, "2.5 %"] <= psi[1] & coefs[2, "97.5 %"] >= psi[1])
  ci_width <- as.numeric(coefs[2, "97.5 %"] - coefs[2, "2.5 %"])
  
  
  ### Placeholder vectors
  #----------------------
  chemical_estimates <- numeric(7)
  ci_lowers <- numeric(7)
  ci_uppers <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  parameters <- c("Joint_Increase_All_Exposures", "X1","X3", "X2", "X7", "X6", "X5", "X4")
  
  
  ### Calculate statistics for each chemical
  #-----------------------------------------
  for(j in seq(1,7,1)){
    coef_w <- resultwqs$final_weights[parameters[j + 1],]
    chemical_estimates[j] <- as.numeric(coef_w[2])
    
    ci <- as.numeric(coef_w[3:4])
    ci_lowers[j] <- ci[1]
    ci_uppers[j] <- ci[2]
    
    ci_widths[j] <- ci[2] - ci[1]
    empirical_powers[j] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[j] <- (ci[1] <= psi[j+1]) & (ci[2] >= psi[j+1])
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
    Method = paste0("WQS (q=",q, ") regression (rh=", rh, " , bs=", bs, ")"),
    Size = nrow(obs),
    Power = empirical_powers_all,
    CI_coverage = CI_coverages_all,
    CI_width = ci_widths_all,
    Correlation = rho
  )
  return(coef_WQS)
}


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Set seed
#-----------
set.seed(123)


### Define parameters
#--------------------
mu <- c(2.1223 , -1.9563 , -1.960, -0.5490 , -1.3725, 0.13153,1.1773 )
correlations <- c(0,0.8)


### Store results
#----------------
scenarioa <- data.frame()


for (n in sample_sizes) {
  for (rho in correlations) {
    
    
    ### Define covariance matrix
    #---------------------------
    Sigma <- diag(p)
    Sigma[1:7, 1:7] <- rho
    diag(Sigma) <- 1
    
    
    ### Repeat this procedure x times
    #--------------------------------
    for (i in 1:n_sim) {
      
      print(i)
      ### Simulate data
      #----------------
      X <- rmvt(n, df = 2, sigma = Sigma)
      X <- sweep(X , 2 , mu, "+")
      X <- scale(X)
      colnames(X) <- c("X1","X3", "X2", "X7", "X6", "X5", "X4")
      X <- as.data.frame(X)
      
      
      ### Sample the error terms using a variance of 0.33 from classical shrinkage methods
      #-----------------------------------------------------------------------------------
      eps <- rnorm(n , mean = 0, sd = sqrt(0.8466327))
    
    
      ### Use a WQS inspired model to simulate the outcome
      #---------------------------------------------------
      WQS <- 0.186305814*X$X1 + 0.051494073*X$X2 + 0.421321419  * X$X3 + 0.020734543 * X$X4 + 0.232626663  * X$X5 + 0.002444214  * X$X6 + 0.085073273  * X$X7
    
    
      ### Combine as a linear model functional form
      #--------------------------------------------
      Y <- -1.471401 + psi[1] * WQS + eps
      obs <- data.frame(Y , X)
      
    
    ### Merge results from all models
    #--------------------------------
    coef_WQS_rh <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                alpha = 0.05, sim_id = i, psi, q = NULL, rho)
      
    coef_WQS_rh_4 <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                  alpha = 0.05, sim_id = i, psi, q = 4, rho)
      
    coef_WQS_rh_10 <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                   alpha = 0.05, sim_id = i, psi, q = 10, rho)
    scenarioa <- rbind(scenarioa , coef_WQS_rh,coef_WQS_rh_4,coef_WQS_rh_10 )
    }
  }
}


### Switch names
#---------------
scenarioa <- scenarioa %>%
  mutate(
    Quantile = case_when(
      Method == "WQS (q=) regression (rh=100 , bs=20)" ~ "Continuous",
      Method == "WQS (q=4) regression (rh=100 , bs=20)" ~ "Quartile",
      Method == "WQS (q=10) regression (rh=100 , bs=20)" ~ "Decile"
    )
  )
scenarioa$Correlation <-  case_when(
  scenarioa$Correlation == 0 ~ "No correlation",
  scenarioa$Correlation == 0.8 ~"High correlation"
    )


### Select coefficients of interest
#----------------------------------
coefwqs <- subset(scenarioa, parameter == "Joint_Increase_All_Exposures")


### Force the order of Quantile levels
#-------------------------------------
coefwqs$Quantile <- factor(coefwqs$Quantile, levels = c("Quartile", "Decile", "Continuous"))
coefwqs$Correlation <- factor(coefwqs$Correlation, levels = c("No correlation",
                                                              "High correlation"))

### Calculate the power
#----------------------
power_summary <- coefwqs %>%
  group_by(Method, Size,Quantile,Correlation) %>%
  summarise(
    Positives = sum(Power == 1, na.rm = TRUE),  # Count of significant results
    n = n(),                                    # Total number of simulations
    Power = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )


### Plot power of the results
#----------------------------
p3<-power_summary %>%
  ggplot(aes(x = Quantile, y = Power, fill = Quantile, color = Quantile)) +  # Add color here
  geom_bar(stat = "identity", alpha = 0.5, position = "dodge") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  facet_grid(Correlation~ " ")+
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 16),
    axis.title.y = element_text(size = 20)
  ) +
  ylab("Empirical power of the joint effect") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ', color = 'Repeated holdout WQS regression with exposure type: ')


### Store the results for each exposure
#--------------------------------------
coefPFAS <- subset(scenarioa, parameter == "X1" | parameter == "X2" |
                 parameter == "X3" |parameter == "X4" |
                parameter == "X5" |parameter == "X6"|parameter == "X7")


### Force the order of Quantile levels
#-------------------------------------
coefPFAS$Quantile <- factor(coefPFAS$Quantile, levels = c("Quartile", "Decile", "Continuous"))
coefPFAS$Correlation <- factor(coefPFAS$Correlation, levels = c("No correlation", "High correlation"))


### Calculate the bias
#---------------------
coefPFAS <- coefPFAS %>%
  mutate(
    # Match the chemical index in the psi vector
    psi_index = case_when(
      parameter == "X1" ~ 2,
      parameter == "X3" ~ 3,
      parameter == "X2" ~ 4,
      parameter == "X7" ~ 5,
      parameter == "X6" ~ 6,
      parameter == "X5" ~ 7,
      parameter == "X4" ~ 8,
    ),
    # Calculate the true effect based on method group
    TrueEffect = psi[psi_index]
    ,
    # Calculate the absolute bias and relative bias
    AbsoluteBias = Estimate - TrueEffect,
    RelativeBias = (Estimate - TrueEffect) / abs(TrueEffect)
  )


### Change name labels
#---------------------
name_map <- c("X1" = "X1 (18%)", "X2" = "X2 (5%)", 
              "X3" = "X3 (42%)", "X4" = "X4 (2%)",
              "X5"="X5 (23%)","X6"="X6 (~0%)", "X7" = "X7 (9%)")
coefPFAS$parameter<-name_map[coefPFAS$parameter]
coefPFAS$parameter<-factor(coefPFAS$parameter , name_map)


### Plot the weights
#-------------------
p4 <- ggplot(coefPFAS, aes(x = Quantile, y = AbsoluteBias, fill = Quantile)) +
  geom_boxplot(alpha=0.5, show.legend = FALSE) +
  geom_jitter(aes(color = Quantile), width = 0.2, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  ylim(-0.4,0.4) + 
  facet_grid(Correlation~ parameter)+
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 20)
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
  ylab("Absolute bias") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ',
       color = 'Repeated holdout WQS regression with exposure type: ') +
  geom_hline(yintercept = 0, col = "red")


########################### Exposure scenario (b) ##############################
#-------------------------------------------------------------------------------
# Exposures from a multivariate normal distribution with added exposure-driven out-
# liers by multiplicative inflation on 5 randomly selected values with either no corre-
#  lation or high correlation among the exposures.
#-------------------------------------------------------------------------------


### Set seed
#-----------
set.seed(123)


### Store results
#----------------
resultsoutlier <- data.frame()


for (n in sample_sizes) {
  for (rho in correlations) {
    
    
    ### Define covariance matrix
    #---------------------------
    Sigma <- diag(p)
    Sigma[1:7, 1:7] <- rho
    diag(Sigma) <- 1
    
    
    ### Repeat this procedure x times
    #--------------------------------
    for (i in 1:n_sim) {
      print(i)
      
      ### Simulate data
      #----------------
      X <- mvrnorm(n, mu = mu, Sigma = Sigma)
      
      
      ### Add outliers
      #---------------
      index <- sample(1:n , 5)
      for(a in index){
        X[a,] <- X[a,] * (2.5 + runif(7,min = 0,max = 1))
      }
      
      
      ### Scale and name
      #-----------------
      X <- scale(X)
      colnames(X) <- c("X1","X3", "X2", "X7", "X6", "X5", "X4")
      X <- as.data.frame(X)
      
      
      ### Sample the error terms using a variance of 0.33 from classical shrinkage methods
      #-----------------------------------------------------------------------------------
      eps <- rnorm(n , mean = 0, sd = sqrt(0.8466327))
    
    
      ### Use a WQS inspired model to simulate the outcome
      #---------------------------------------------------
      WQS <- 0.186305814*X$X1 + 0.051494073*X$X2 + 0.421321419  * X$X3 + 
        0.020734543 * X$X4 + 0.232626663  * X$X5 + 0.002444214  * X$X6 + 0.085073273  * X$X7
    
    
      ### Combine as a linear model functional form
      #--------------------------------------------
      Y <- -1.471401 + psi[1] * WQS + eps
      
      
      ### Store the simulated data
      #---------------------------
      obs <- data.frame(Y , X)
      
    
      ### Fit the three different models 
      #---------------------------------
      coef_WQS_rh <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                alpha = 0.05, sim_id = i, psi, q = NULL, rho)
      
      coef_WQS_rh_4 <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                  alpha = 0.05, sim_id = i, psi, q = 4, rho)
      
      coef_WQS_rh_10 <- run_WQS_rh(obs, rh = 100, bs = 20, signal = "one" , 
                                   alpha = 0.05, sim_id = i, psi, q = 10, rho)
  
      
      ### Merge results from all models
      #--------------------------------
       resultsoutlier <- rbind(resultsoutlier,coef_WQS_rh,coef_WQS_rh_4,coef_WQS_rh_10 )
    }
  }
}


### Switch names
#---------------
resultsoutlier <- resultsoutlier %>%
  mutate(
    Quantile = case_when(
      Method == "WQS (q=) regression (rh=100 , bs=20)" ~ "Continuous",
      Method == "WQS (q=4) regression (rh=100 , bs=20)" ~ "Quartile",
      Method == "WQS (q=10) regression (rh=100 , bs=20)" ~ "Decile"
    )
  )
resultsoutlier$Correlation <-  case_when(
      resultsoutlier$Correlation == 0 ~ "No correlation",
      resultsoutlier$Correlation == 0.8 ~"High correlation"
    )


### Store joint effect
#---------------------
coefwqsoutlier <- subset(resultsoutlier, parameter == "Joint_Increase_All_Exposures")



### Force the order of Quantile levels
#-------------------------------------
coefwqsoutlier$Quantile <- factor(coefwqsoutlier$Quantile, levels = c("Quartile", "Decile", "Continuous"))
coefwqsoutlier$Correlation <- factor(coefwqsoutlier$Correlation, levels = c("No correlation",
                                                              "High correlation"))

### Calculate the power
#----------------------
power_summary_outlier <- coefwqsoutlier %>%
  group_by(Method, Size,Quantile,Correlation) %>%
  summarise(
    Positives = sum(Power == 1, na.rm = TRUE),  # Count of significant results
    n = n(),                                    # Total number of simulations
    Power = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )


### Plot the empirical power
#---------------------------
p5<-power_summary_outlier %>%
  ggplot(aes(x = Quantile, y = Power, fill = Quantile, color = Quantile)) +  # Add color here
  geom_bar(stat = "identity", alpha = 0.5, position = "dodge") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  facet_grid(Correlation~ " ")+
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 16),
    axis.title.y = element_text(size = 20)
  ) +
  ylab("Empirical power of the joint effect") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ', color = 'Repeated holdout WQS regression with exposure type: ')



### Select all PFAS components
#-----------------------------
coefPFASoutlier <- subset(resultsoutlier, parameter == "X1" | parameter == "X2" |
                 parameter == "X3" |parameter == "X4" |
                parameter == "X5" |parameter == "X6"|parameter == "X7")


### Force the order of Quantile levels
#-------------------------------------
coefPFASoutlier$Quantile <- factor(coefPFASoutlier$Quantile, 
                            levels = c("Quartile" , "Decile" , "Continuous"))
coefPFASoutlier$Correlation <- factor(coefPFASoutlier$Correlation, 
                               levels = c("No correlation" , "High correlation"))


### Calculate the bias
#---------------------
coefPFASoutlier <- coefPFASoutlier %>%
  mutate(
    # Match the chemical index in the psi vector
    psi_index = case_when(
      parameter == "X1" ~ 2,
      parameter == "X3" ~ 3,
      parameter == "X2" ~ 4,
      parameter == "X7" ~ 5,
      parameter == "X6" ~ 6,
      parameter == "X5" ~ 7,
      parameter == "X4" ~ 8,
    ),
    # Calculate the true effect based on method group
    TrueEffect = psi[psi_index]
    ,
    # Calculate the absolute bias and relative bias
    AbsoluteBias = Estimate - TrueEffect,
    RelativeBias = (Estimate - TrueEffect) / abs(TrueEffect)
  )


### Change name labels
#---------------------
name_map <- c("X1" = "X1 (18%)", "X2" = "X2 (5%)", 
              "X3" = "X3 (42%)", "X4" = "X4 (2%)",
              "X5"="X5 (23%)","X6"="X6 (~0%)", "X7" = "X7 (9%)")
coefPFASoutlier$parameter<-name_map[coefPFASoutlier$parameter]
coefPFASoutlier$parameter<-factor(coefPFASoutlier$parameter , name_map)



### Plot the bias of the weights
#-------------------------------
p6 <- ggplot(coefPFASoutlier, aes(x = Quantile, y = AbsoluteBias, fill = Quantile)) +
  geom_boxplot(alpha=0.5, show.legend = FALSE) +
  geom_jitter(aes(color = Quantile), width = 0.2, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  facet_grid(Correlation~ parameter)+
  theme_ipsum() +
  ylim(-0.4,0.4)+
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 14),
    axis.title.y = element_text(size = 20)
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
  ylab("Absolute bias") +
  labs(fill = 'Repeated holdout WQS regression with exposure type: ',
       color = 'Repeated holdout WQS regression with exposure type: ') +
  geom_hline(yintercept = 0, col = "red")



################################################################################
############################ Combine all plots #################################
################################################################################


### Row titles
#-------------
row1_title <- wrap_elements(full = grid::textGrob("Exposure scenario (a): Multivariate t-distribution (df=2) ", gp = grid::gpar(fontsize = 16, fontface = "bold")))
row2_title <- wrap_elements(full = grid::textGrob("Exposure scenario (b): Multivariate normal distribution with 5 exposure-outliers", gp = grid::gpar(fontsize = 16, fontface = "bold")))
row3_title <- wrap_elements(full = grid::textGrob("Exposure scenario (c):  Resampling with replacement", gp = grid::gpar(fontsize = 16, fontface = "bold")))


### Combine plots with row titles
#--------------------------------
final_plot <- (
  row1_title / (p4 + p3 + plot_layout(widths = c(0.93, 0.07))) /
  row2_title / (p6 + p5 + plot_layout(widths = c(0.93, 0.07))) /
  row3_title / (p2 + p1 + plot_layout(widths = c(0.93, 0.07))) +
  plot_layout(heights = c(0.005, 0.37,0.005, 0.37,0.006, 0.26), guides = "collect")
) & theme(legend.position = "bottom")

