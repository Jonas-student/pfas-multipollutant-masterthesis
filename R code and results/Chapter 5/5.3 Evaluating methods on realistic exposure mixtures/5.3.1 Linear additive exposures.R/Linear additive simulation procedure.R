################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Linear additive exposures
# 24/06/2025
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
library(doParallel)
library(parallel)
library(randomForest)
library(bkmr)
library(BAS)
library(glmnet)
library(nimble)

### Load data for a realistic simulation study
#---------------------------------------------
df <- read_excel("~/Dataset_HBM_3M.xlsx")


### Make gender a dummy variable (0/1)
#-------------------------------------
df$GESL <- (df$GESL-1)


### Create data frame
#--------------------
Data <- data.frame("lbpfos" = log(df$lbpfos_imp) , "lbpfhxs" = log(df$lbpfhxs_imp), "lbpfoa"=log(df$lbpfoa_imp),
                   "pfba" = log(df$pfba_imp), "pfda" = log(df$pfda_imp) , 
                   "pfna" = log(df$pfna_imp) , "pfos" = log(df$pfos_imp),
                   "age"=df$OB_Leeftijd, "highest education"=df$ISCED_HH , "gender"=df$GESL)


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


### Sample size and number of variables
#--------------------------------------
n <- nrow(data_matrix)
p <- ncol(data_matrix)


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Write a parallel simulation function
#---------------------------------------
estimates <-  function(m, seed){
  
  
  ### Define the effect size
  #-------------------------
  psi <- c(-0.3 , 0.186305814, 0.421321419,  0.051494073, 0.085073273, 0.002444214,
           0.232626663, 0.020734543)
  alpha <- 0.05
  
  
  ### Function to simulate data
  #----------------------------
  sim <- function(n_sim){
    
    
    ### Resample indices (bootstrap resampling)
    #------------------------------------------
    boot_indices <- sample(1:300, size = n_sim, replace = TRUE)
    
    
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
    eps <- rnorm(n_sim , mean = 0, sd = sqrt(0.8466327/4))
    
    
    ### Use a WQS inspired model to simulate the outcome
    #---------------------------------------------------
    WQS <- 0.186305814*X$lbpfos + 0.051494073*X$pfda + 0.421321419* X$pfba + 
      0.020734543*X$pfos + 0.232626663*X$lbpfoa + 0.002444214*X$pfna + 0.085073273*X$lbpfhxs
    
    
    ### Combine as a linear model functional form
    #--------------------------------------------
    Y <- -1.471401 + psi[1] * WQS - 0.016248 * X$age + 0.168078 * X$gender + 
      0.103059 * X$Educ2 + 0.023381 * X$Educ3 + eps
    
    return(data.frame(Y , X))
  }
  
  
  ### Parellise for loop
  #---------------------
  out <- foreach(i=1:m, .combine=rbind, .packages = c("gWQS", "MASS","qgcomp","randomForest","bkmr","BAS","glmnet","nimble","bkmr"))%dopar%{
    
    print(i)
    
    ### Set the sample size
    #----------------------
    n_sim <- 300
    q <- NULL
    
    
    ### Set seed
    #-----------
    set.seed(seed[i])
    
    
    ### Simulated observed data
    #--------------------------
    obs <- sim(n_sim)
    
    
    ### Quantile g-computation package
    #---------------------------------
    source("~/gcomp_function_store_results.R")
    coef_gcomp <- run_qgcomp(obs, bs = 200, alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### Quantile g-computation linear model
    #---------------------------------------
    source("~/lm_gcomp_function_store_results.R")
    coef_lm_gcomp <- run_lm_gcomp_boot(obs, exposure_shift = 1, n_boot = 200, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Quantile g-computation random forest
    #---------------------------------------
    source("~/rf_gcomp_function_store_results.R")
    coef_rf <- run_rf_gcomp_boot(obs, exposure_shift = 1, n_boot = 200, alpha = 0.05, ntree = 500, sim_id = seed[i], psi = psi)
    
    
    ### BAS
    #------
    source("~/BAS_function_store_results.R")
    coef_BAS <- run_BAS(obs, alpha = 0.05, sim_id = seed[i], psi = psi)
    
    
    ### Linear model (OLS regression)
    #-------------------------------
    source("~/OLS_function_store_results.R")
    coef_lm <- run_linear_model(obs, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Linear model (OLS regression)
    #-------------------------------
    source("~/singleOLS_function_store_results.R")
    coef_singlelm <- run_singlelinear_model(obs, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Ridge regression
    #-------------------
    source("~/Ridge_function_store_results.R")
    coef_Ridge <- run_Ridge(obs, n_boot = 200, alpha = 0.05, sim_id = seed[i], psi = psi)
    
    
    ### Lasso regression
    #-------------------
    source("~/Lasso_function_store_results.R")
    coef_Lasso <- run_Lasso(obs, n_boot = 200, alpha = 0.05, sim_id = seed[i], psi = psi)
    
    
    
    ### ENET regression
    #------------------
    source("~/ENET_function_store_results.R")
    coef_ENET <- run_enet(obs, n_boot = 200, alpha = 0.05, sim_id = seed[i], psi = psi)
    
    
    ### WQS bootstrap 
    #----------------
    source("~/WQS_bs_function_store_results.R")
    coef_bs <- run_WQS_bs(obs, bs = 100, signal = "one" , alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### WQS abst bootstrap 
    #---------------------
    coef_bs_abst <- run_WQS_bs(obs, bs = 100, signal = "abst" , alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### WQS random subset 
    #--------------------
    source("~/WQS_rs_function_store_results.R")
    coef_rs <- run_WQS_rs(obs, bs = 100, signal = "one" , alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### WQS repeated holdout
    #-----------------------
    source("~/WQS_rh_function_store_results.R")
    coef_WQS_rh <- run_WQS_rh(obs, rh = 100, bs = 10, signal = "one" , alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### WQS repeated holdout apply the absolute value of the t-statistic
    #-------------------------------------------------------------------
    coef_WQS_rh_abst <- run_WQS_rh(obs, rh = 100, bs = 10, signal = "abst" , alpha = 0.05, sim_id = seed[i], psi, q = NULL)
    
    
    ### Bayesian Lasso
    #-----------------
    source("~/BayesianLasso_function_store_results.R")
    coef_BayLasso <- run_BLasso(obs, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Bayesian spike & slab
    #------------------------
    source("~/SpikeSlab_function_store_results.R")
    coef_Spikeslab <- run_SpikeSlab(obs, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Bayesian Horseshoe
    #---------------------
    source("~/Horseshoe_function_store_results.R")
    coef_Horseshoe <- run_Horseshoe(obs, alpha = 0.05, sim_id = seed[i], psi)
    
    
    ### Bayesian kernel machine regression
    #-------------------------------------
    source("~/bkmr_function_store_results.R")
    coef_bkmr <- run_bkmr(obs, sim_id = seed[i], psi)
    
    
    ### Merge results from all models
    #--------------------------------
    rbind( coef_bs, coef_bs_abst,coef_rs, coef_WQS_rh,
           coef_WQS_rh_abst , coef_gcomp, coef_rf, coef_lm_gcomp, coef_lm, coef_BAS, 
           coef_Ridge, coef_Lasso, coef_ENET, coef_BayLasso, coef_Spikeslab,
           coef_Horseshoe, coef_singlelm, coef_bkmr)
  }
  return(out)
}


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Initialise log file
#----------------------
log_file <- "simulation_log.txt"
cat("Simulation started at", as.character(Sys.time()), "\n", file = log_file, append = FALSE)


### Define all seeds
#-------------------
seed <- 85721023 + 0:99


### Check for existing checkpoint
#--------------------------------
checkpoint_file <- "checkpoint_results.rds"
if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  start_batch <- checkpoint$last_batch + 1
  results <- checkpoint$results
  cat("Resuming from batch", start_batch, "\n", file = log_file, append = TRUE)
} else {
  start_batch <- 1
  results <- list()
  cat("Starting new simulation\n", file = log_file, append = TRUE)
}


### Define total simulations
#---------------------------
total_batches <- 25
iterations_per_batch <- 4
total_iterations <- total_batches * iterations_per_batch


### Run simulation with timing
#-----------------------------
timing <- system.time({
  
  for (batch in start_batch:total_batches) {
    
    ### Intermediate check
    #---------------------
    print(batch)
    
    
    ### Define random seed
    #---------------------
    seed_subset <- seed[(4 * (batch - 1) + 1):(4 * batch)]
    
    
    ### Register cluster
    #-------------------
    cluster <- makeCluster(4)
    registerDoParallel(cluster)
    
    
    ### Export necessary variables and functions to the cluster
    #----------------------------------------------------------
    clusterExport(cluster, varlist = c( "Data", "df", "n", "p", "data_matrix","PFAS_clean"), envir = environment())
    
    
    ### Lapse total time
    #-------------------
    system.time({
      
    batch_results <- tryCatch({
        estimates(iterations_per_batch , seed_subset)
      }, 
      error = function(e) {
        NULL
      })
    
    })
    
    
    ### Stop cluster
    #---------------
    stopCluster(cluster)
    closeAllConnections()
    

    ### Store batch results
    #----------------------
    if (!is.null(batch_results)) {
      results[[as.character(batch)]] <- batch_results
      saveRDS(batch_results, file = paste0("batch_", batch, ".rds"))
      cat("Batch", batch, "results saved to batch_", batch, ".rds\n", file = log_file, append = TRUE)
    }
    
    
    ### Save checkpoint
    #------------------
    saveRDS(list(last_batch = batch, results = results), file = checkpoint_file)
    cat("Checkpoint saved at batch", batch, "(iteration", batch * iterations_per_batch, ")\n", file = log_file, append = TRUE)
    
  }
})



