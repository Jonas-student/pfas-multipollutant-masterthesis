################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Histogram with Z-score outliers highlighted distributions
# Appendix C
# 10/07/2025
################################################################################


### Load libraries
#-----------------
library(ggplot2)
library(MASS)


### Define covariance matrix
#---------------------------
Sigma <- diag(p)
Sigma[1:7, 1:7] <- 0.5
diag(Sigma) <- 1


### Set seed for reproducibility
#-------------------------------
set.seed(1234)


### Generate exposures from a multivariate normal distribution
#-------------------------------------------------------------
X <- mvrnorm(n, mu = mu, Sigma = Sigma)


### Add outliers
#---------------
index <- sample(1:n , 5)
for(i in index){
  X[i,] <- X[i,] * (2.5 + runif(7,min = 0,max = 1))
}


### Scale and name
#-----------------
X <- scale(X)
colnames(X) <- c("lbpfos", "pfba", "pfda", "lbpfhxs", "pfna", "lbpfoa", "pfos")
X <- as.data.frame(X)


### Compute Z-scores
#-------------------
z_scores <- scale(X$lbpfos)


### Create a data frame with an outlier flag
#-------------------------------------------
df_hist <- data.frame(x = X$lbpfos, z = z_scores)
df_hist$outlier <- ifelse(abs(df_hist$z) > 3, "Outlier", "Normal")


### Plot histogram
#-----------------
ggplot(df_hist, aes(x = x, fill = outlier)) +
  geom_histogram(bins = 30, alpha = 0.75, color = "black") +
  scale_fill_manual(values = c("Normal" = "gray", "Outlier" = "red")) +
  labs(title = "Histogram with Z-score outliers highlighted",
       x = "Exposure Value", y = "Count",fill = "Exposure type:") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_blank(),        # Remove x-axis labels
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) 
