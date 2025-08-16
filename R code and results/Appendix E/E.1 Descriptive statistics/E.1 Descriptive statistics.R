################################################################################
# Case study - Jonas Meijerink
# Descriptive statistics
# 27/07/2025
################################################################################


rm(list = ls())


### Load libraries 
#-----------------
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmnet)


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
Data <- data.frame("PFOS (branched)" = log(df$bpfos_imp) , "PFHXS (total)" = log(df$lbpfhxs_imp), 
                   "PFOA (total)"=log(df$lbpfoa_imp), "PFBA" = log(df$pfba_imp), 
                   "PFDA" = log(df$pfda_imp) , "PFNA" = log(df$pfna_imp) , 
                   "PFOS" = log(df$pfos_imp), "age2"=df$Age2, "age3"=df$Age3,
                   "gender"=df$GESL, "income2"=df$Income2,
                   "income3"=df$Income3 , "bmi2"=df$BMI2, "bmi3"=df$BMI3,
                   "smoke"=df$ALLSMOKE, "birthweight"=df$GEBOORTEGEW_LICHT, 
                   "lbpfos_I"=(df$lbpfos==-3),"lbpfhxs_I"=(df$lbpfhxs==-3),"lbpfoa_I"=(df$lbpfoa==-3),
                   "pfba_I"=(df$pfba==-3),"pfda_I"=(df$pfda==-3),"pfna_I"=(df$pfna==-3),
                   "pfos_I"=(df$pfos==-3))
Data[1:7] <- scale(Data[1:7])


### Define the outcome
#---------------------
Y <- log(df$cd4plus_imp/df$cd8plus_imp)
Y2 <- log(df$leuko)


### Delete missing values
#------------------------
obs <- data.frame(Y , Data)
obs <- na.omit(obs)

obs2 <- data.frame(Y2 , Data)
obs2 <- na.omit(obs2)


################################################################################
#-------------------------- Descriptive statistics ----------------------------#
################################################################################


### Descriptive statistics summary table (Table E.1)
#---------------------------------------------------
summary(df$pfos)
summary(df$bpfos_imp)
summary(df$lbpfhxs)
summary(df$lbpfoa)
summary(df$pfba)
summary(df$pfda)


### Visualize exposure distribution
#----------------------------------
plots <- list()
names <- c("PFOS (branched)", "PFHXS (total)", 
           "PFOA (total)", "PFBA","PFDA", "PFNA", "PFOS")

for (i in 2:8) {
  x <- obs[, i]
  below_loq <- obs[, i + 16]
  
  df_plot <- data.frame(x = x, below_loq = below_loq)
  
  df_plot$below_loq <- factor(
    df_plot$below_loq,
    levels = c(TRUE, FALSE),
    labels = c("< LOQ", "≥ LOQ")
  )
  
  n_below <- sum(df_plot$below_loq == "< LOQ", na.rm = TRUE)
  n_above <- sum(df_plot$below_loq == "≥ LOQ", na.rm = TRUE)
  
  p <- ggplot(df_plot, aes(x = x)) +
    # Histogram with light gray fill, semi-transparent
    geom_histogram(aes(y = ..density..), bins = 20, fill = "grey90", color = "white", alpha = 1) +
    
    # Density for < LOQ
    geom_density(data = subset(df_plot, below_loq == "< LOQ"),
                 aes(color = below_loq),
                 size = 0.9, adjust = 1.5, linetype = "dashed") +
    
    # Density for ≥ LOQ
    geom_density(data = subset(df_plot, below_loq == "≥ LOQ"),
                 aes(color = below_loq),
                 size = 0.9, adjust = 1.5) +
    
    scale_color_manual(
      values = c("< LOQ" = "tomato", "≥ LOQ" = "skyblue"),
      labels = c(paste0("< (n=", n_below, ")"), paste0("≥ (n=", n_above, ")")),
      name = "LOQ status: ",
      drop = FALSE
    ) +
    
    labs(
      x = paste("log(", names[i-1], ")", sep = ""),
      y = "Density"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.ticks = element_line(color = "grey80", size = 0.3),
      axis.text = element_text(size = 14)
    )
  
  plots[[i - 1]] <- p
}

# Combine plots
wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Histogram of log-transformed chemical distributions with density curve by LOQ",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )


### Distribution of the outcome
#------------------------------
p1 <- ggplot(obs,aes(x=Y)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "grey90", color = "white", alpha = 1) +
  geom_density(aes(y = ..density..),
               size = 0.9, adjust = 1.5) +
  labs(
    x = "log(CD4+/CD8+ T-cell ratio)",
    y = "Density"
  )+
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.ticks = element_line(color = "grey80", size = 0.3),
    axis.text = element_text(size = 14)
  )
p2 <- ggplot(obs2,aes(x=Y2)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "grey90", color = "white", alpha = 1) +
  geom_density(aes(y = ..density..),
               size = 0.9, adjust = 1.5) +
  labs(
    x = "log(Leukocyte count)",
    y = "Density"
  )+
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.ticks = element_line(color = "grey80", size = 0.3),
    axis.text = element_text(size = 14)
  )
p1+p2



################################################################################
#                                                                              #
#----------------------------- Leukocyte count --------------------------------#
#                                                                              #
################################################################################


### Center and scale exposures
#-----------------------------
Data <- data.frame("bpfos" = log(df$bpfos_imp) , "lbpfhxs" = log(df$lbpfhxs_imp), 
                   "lbpfoa"=log(df$lbpfoa_imp), "pfba" = log(df$pfba_imp), 
                   "pfda" = log(df$pfda_imp) , "pfna" = log(df$pfna_imp) , 
                   "pfos" = log(df$pfos_imp), "age2"=df$Age2, "age3"=df$Age3,
                   "gender"=df$GESL, "income2"=df$Income2,
                   "income3"=df$Income3 , "bmi2"=df$BMI2, "bmi3"=df$BMI3,
                   "smoke"=df$ALLSMOKE, "birthweight"=df$GEBOORTEGEW_LICHT)
Data[1:7] <- scale(Data[1:7])


### Define the outcome
#---------------------
Y <- log(df$leuko)


### Delete missing values
#------------------------
obs <- data.frame(Y , Data)
obs <- na.omit(obs)


### Data frame with log-transformed response and predictors
#----------------------------------------------------------
df <- data.frame(y = obs$Y, 
                 x1 = obs$bpfos, x2 = obs$lbpfhxs, x3 = obs$lbpfoa, 
                 x4 = obs$pfba, x5 = obs$pfda, x6 = obs$pfna, x7 = obs$pfos)


### List of predictor names
#--------------------------
predictors <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7")
predictor_labels <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", "pfos")


### Generate all two-way combination plots
#-----------------------------------------
for (i in 1:(length(predictors)-1)) {
  for (j in (i+1):length(predictors)) {
    p <- ggplot(df, aes_string(x = predictors[i], y = "y")) +
      geom_point(alpha = 0.5) +
      facet_wrap(as.formula(paste("~ cut(", predictors[j], ", 3)")), ncol = 1) +
      geom_smooth( se = TRUE, color = "blue") +
      labs(title = paste("log(Y) vs", predictor_labels[i], "by", predictor_labels[j]),
           x = predictor_labels[i], y = "log(Y)") +
      theme_minimal()
    print(p)
  }
}


### Full model with pair-wise interactions
#-----------------------------------------
full_model <- Y ~ bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos +
  bpfos:lbpfhxs + bpfos:lbpfoa + bpfos:pfba + bpfos:pfda + bpfos:pfna + bpfos:pfos +
  lbpfhxs:lbpfoa + lbpfhxs:pfba + lbpfhxs:pfda + lbpfhxs:pfna + lbpfhxs:pfos +
  lbpfoa:pfba + lbpfoa:pfda + lbpfoa:pfna + lbpfoa:pfos + 
  pfba:pfda + pfba:pfna + pfba:pfos +
  pfda:pfna + pfda:pfos +
  pfna:pfos +
  age2 + age3 + gender + income2 + income3 + bmi2 + bmi3 + smoke + birthweight


### Prepare data
#---------------
X <- model.matrix(full_model, data = obs)[, -1]
Y <- obs$Y
penalty_factors <- c(rep(0, 7), rep(0, 9), rep(1 , 21))
lambda <- exp(seq(log(1e-6), log(1e3), length.out = 10000))


### Fit Lasso coefficients path
#------------------------------
lasso_fit <- glmnet(X, Y, alpha = 1, penalty.factor = penalty_factors,
                    lambda = lambda)
coef_mat <- as.matrix(lasso_fit$beta)[c(17:37), ]  # exclude intercept
coef_mat <- coef_mat[, ncol(coef_mat):1]
matplot(log10(lambda), t(coef_mat), type = "l", lty = 1, col = 1:7, 
        xlab = "log10(Lambda)", ylab = "Coefficient", main = "Lasso coefficients")
abline(v = -2, lty = 2)
abline(v = -4, lty = 2)


### Bootstrap inclusion probabilities
#------------------------------------
set.seed(485721023)
n_boot <- 2000
stability <- matrix(0, nrow = 21, ncol = n_boot)  
for (b in 1:n_boot) {
  boot_idx <- sample(1:nrow(obs), replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X, Y, alpha = 0, penalty.factor = penalty_factors,
                      type.measure = "mse", nfolds = 5, lambda = 
                        exp(seq(log(1e-4), log(1e-2), length.out = 1000)))
  fit_boot <- glmnet(X_boot, Y_boot, alpha = 1, penalty.factor = penalty_factors,
                     lambda = cv_fit$lambda.1se)
  coefs <- as.matrix(coef(fit_boot))[-1, ]  
  stability[, b] <- coefs[17:37]!=0  
}
stability_scores <- rowMeans(stability)
names(stability_scores) <- colnames(X)[17:37]  
print(stability_scores)  




################################################################################
#                                                                              #
#----------------------------- CD4+/CD8+ ratio --------------------------------#
#                                                                              #
################################################################################


### Define the outcome
#---------------------
Y <- log(df$cd4plus_imp/df$cd8plus_imp)


### Delete missing values
#------------------------
obs <- data.frame(Y , Data)
obs <- na.omit(obs)


### Data frame with log-transformed response and predictors
#----------------------------------------------------------
df <- data.frame(y = obs$Y, 
                 x1 = obs$bpfos, x2 = obs$lbpfhxs, x3 = obs$lbpfoa, 
                 x4 = obs$pfba, x5 = obs$pfda, x6 = obs$pfna, x7 = obs$pfos)


### List of predictor names
#--------------------------
predictors <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7")
predictor_labels <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", "pfos")


### Generate all two-way combination plots
#-----------------------------------------
for (i in 1:(length(predictors)-1)) {
  for (j in (i+1):length(predictors)) {
    p <- ggplot(df, aes_string(x = predictors[i], y = "y")) +
      geom_point(alpha = 0.5) +
      facet_wrap(as.formula(paste("~ cut(", predictors[j], ", 3)")), ncol = 1) +
      geom_smooth( se = TRUE, color = "blue") +
      labs(title = paste("log(Y) vs", predictor_labels[i], "by", predictor_labels[j]),
           x = predictor_labels[i], y = "log(Y)") +
      theme_minimal()
    print(p)
  }
}


### Full model with pair-wise interactions
#-----------------------------------------
full_model <- Y ~ bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos +
  bpfos:lbpfhxs + bpfos:lbpfoa + bpfos:pfba + bpfos:pfda + bpfos:pfna + bpfos:pfos +
  lbpfhxs:lbpfoa + lbpfhxs:pfba + lbpfhxs:pfda + lbpfhxs:pfna + lbpfhxs:pfos +
  lbpfoa:pfba + lbpfoa:pfda + lbpfoa:pfna + lbpfoa:pfos + 
  pfba:pfda + pfba:pfna + pfba:pfos +
  pfda:pfna + pfda:pfos +
  pfna:pfos +
  age2 + age3 + gender + income2 + income3 + bmi2 + bmi3 + smoke + birthweight


### Prepare data
#---------------
X <- model.matrix(full_model, data = obs)[, -1]
Y <- obs$Y
penalty_factors <- c(rep(0, 7), rep(0, 9), rep(1 , 21))
lambda <- exp(seq(log(1e-6), log(1e3), length.out = 10000))


### Fit Lasso coefficients path
#------------------------------
lasso_fit <- glmnet(X, Y, alpha = 1, penalty.factor = penalty_factors,
                    lambda = lambda)
coef_mat <- as.matrix(lasso_fit$beta)[c(17:37), ]  # exclude intercept
coef_mat <- coef_mat[, ncol(coef_mat):1]
matplot(log10(lambda), t(coef_mat), type = "l", lty = 1, col = 1:7, 
        xlab = "log10(Lambda)", ylab = "Coefficient", main = "Lasso coefficients")
abline(v = -2, lty = 2)
abline(v = -4, lty = 2)


### Bootstrap inclusion probabilities
#------------------------------------
set.seed(485721023)
n_boot <- 2000
stability <- matrix(0, nrow = 21, ncol = n_boot)  
for (b in 1:n_boot) {
  boot_idx <- sample(1:nrow(obs), replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X, Y, alpha = 0, penalty.factor = penalty_factors,
                      type.measure = "mse", nfolds = 5, lambda = 
                        exp(seq(log(1e-4), log(1e-2), length.out = 1000)))
  fit_boot <- glmnet(X_boot, Y_boot, alpha = 1, penalty.factor = penalty_factors,
                     lambda = cv_fit$lambda.1se)
  coefs <- as.matrix(coef(fit_boot))[-1, ]  
  stability[, b] <- coefs[17:37]!=0  
}
stability_scores <- rowMeans(stability)
names(stability_scores) <- colnames(X)[17:37]  
print(stability_scores)  