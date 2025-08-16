################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Intuition behind joint effect SE
# Appendix B
# 23/06/2025
################################################################################

# A common issue with mixtures is that they are highly correlated. 
# A one-unit change in one exposure, holding the other exposure constant, 
# may not be well supported by the data and depend heavily on model extrapolation.


### Load packages
#----------------
library(MASS)
library(ggplot2)


### Sample from a multivariate normal distribution
#-------------------------------------------------
Sigma <- diag(2)
Sigma[1:2, 1:2] <- 0.9
diag(Sigma) <- 1
X <- mvrnorm(300, mu = rep(0,2), Sigma = Sigma)


### Plot Exposure A vs Exposure B
#--------------------------------
plot(X[,1], X[,2], xlab = "Exposure A", ylab = "Exposure B", alpha=0.5)


### Arrow from a reference point
#-------------------------------
x0 <- mean(X[,1], na.rm = TRUE)-0.5
y0 <- mean(X[,2], na.rm = TRUE)-0.5


### Plot arrows to indicate the effect
#-------------------------------------
arrows(x0, y0, x0 + 1, y0, col = "#0072B2", length = 0.1,, lwd = 3, cex=1.5)      
text(x0 + 1, y0, "Exposure A + 1", pos = 4, col = "#0072B2", lwd = 3, cex=1.5)
arrows(x0, y0, x0, y0 + 1, col = "#009E73", length = 0.1, lwd = 3, cex=1.5)     
text(x0, y0 + 1.1, "Exposure B + 1", pos = 3, col = "#009E73", lwd = 3, cex=1.5)
arrows(x0, y0, x0 + 1, y0 + 1, col = "#D55E00", length = 0.1, lwd = 3, cex=1.5)
text(x0 + 1.09, y0 + 1-0.3, "Joint exposure + 1", pos = 4, col = "#D55E00", lwd = 3, cex=1.5)

