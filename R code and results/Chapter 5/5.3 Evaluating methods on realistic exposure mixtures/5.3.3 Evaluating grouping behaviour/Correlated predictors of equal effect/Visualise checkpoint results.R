################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Linear additive exposures visualisations
# 10/07/2025
# 
# Visualisation for selection properties of with two equal effect among two 
# highly correlated predictors
################################################################################


library(ggpubr)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(dplyr)


### Load previous results
#-------------------------
output <- do.call(rbind, checkpoint_results$results)


### Filter on WQS estimate
#-------------------------
coefwqs <- subset(output, parameter == "wqs"|parameter == "WQS"| parameter == "Joint_Increase_All_Exposures" | parameter == "Sum betas")
coefwqs$Estimate <- coefwqs$Estimate + 0.3
coefwqs <- subset(coefwqs, Method != "G-computation (package) (bs=200)")


### Define method order
#----------------------
method_order <- c(
  "Single pollutant linear model",
  "Multiple pollutant linear model",
  
  "Ridge regression (bs=200)",
  "Lasso regression (bs=200)",
  "Elastic Net regression (bs=200)",
  "Bayesian Lasso regression",
  "Bayesian spike & slab regression",
  "Bayesian horseshoe regression",
  
  "G-computation (package) (bs=200)",
  "Linear model g-computation (bs=200)",
  "Random forest g-computation (bs=200)",
  "Bayesian kernel machine regression",
  
  "Bayesian model averaging",
  
  "WQS (one) regression (rh=1 , bs=100)",
  "WQS (abst) regression (rh=1 , bs=100)",
  "WQS (one) regression (100 random subsets)",
  "WQS (one) regression (rh=100 , bs=10)",
  "WQS (abst) regression (rh=100 , bs=10)"
)


### Define group order
#---------------------
group_order <- c(
  "OLS methods",
  "Shrinkage methods",
  "G-computation",
  "Model averaging",
  "Weighted index",
  "Semi-parametric"
)


### Generate a color vector of the same length
#---------------------------------------------
method_colors <- viridis(length(method_order), option = "C", alpha = 0.7)
names(method_colors) <- method_order


### Assign groups
#----------------
coefwqs <- coefwqs %>%
  mutate(Method = factor(Method, levels = method_order),  # Correct method order
         Group = case_when(
           Method %in% c("Multiple pollutant linear model","Single pollutant linear model") ~ "OLS methods",
           Method %in% c("Bayesian Lasso regression", "Elastic Net regression (bs=200)",
                         "Lasso regression (bs=200)", "Ridge regression (bs=200)",
                         "Bayesian spike & slab regression", "Bayesian horseshoe regression") ~ "Shrinkage methods",
           Method %in% c("G-computation (package) (bs=200)","Bayesian kernel machine regression",
                         "Linear model g-computation (bs=200)", "Random forest g-computation (bs=200)") ~ "G-computation",
           Method == "Bayesian model averaging" ~ "Model averaging",
           TRUE ~ "Weighted index"
         ),
         Group = factor(Group, levels = group_order))  # Correct group order


################################################################################

############################# PFAS exposures ###################################

################################################################################




psi <- c(-0.3 , 0, 0,  0.5, 0, 0.5, 0, 0)


### Filter on PFAS estimate
#-------------------------
coefPFAS <- subset(output, parameter == "lbpfos" | parameter == "pfba" |
                     parameter == "pfda" |parameter == "lbpfhxs" |
                     parameter == "pfna" |parameter == "lbpfoa"|parameter == "pfos")


### Assign groups
#----------------
coefPFAS <- coefPFAS %>%
  mutate(Method = factor(Method, levels = method_order),  # Correct method order
         Group = case_when(
           Method %in% c("Multiple pollutant linear model","Single pollutant linear model") ~ "OLS methods",
           Method %in% c("Bayesian Lasso regression", "Elastic Net regression (bs=200)",
                         "Lasso regression (bs=200)", "Ridge regression (bs=200)",
                         "Bayesian spike & slab regression", "Bayesian horseshoe regression") ~ "Shrinkage methods",
           Method %in% c("G-computation (package) (bs=200)",
                         "Linear model g-computation (bs=200)", "Random forest g-computation (bs=200)") ~ "G-computation",
           Method == "Bayesian model averaging" ~ "Model averaging",
           Method %in% c("Bayesian kernel machine regression")~"Semi-parametric",
           TRUE ~ "Weighted index"
         ),
         Group = factor(Group, levels = group_order))  # Correct group order


### Compute the bias
#-------------------
# Compute bias depending on the method group
coefPFAS <- coefPFAS %>%
  filter(Method != "G-computation (package) (bs=200)") %>%
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
    
    # Calculate the true effect based on method group
    TrueEffect = ifelse(Group == "Weighted index",
                        psi[psi_index],            # WQS uses individual effect
                        psi[1] * psi[psi_index]    # Other methods use joint effect
    ),
    
    # Calculate the absolute bias and relative bias
    AbsoluteBias = Estimate - TrueEffect,
    RelativeBias = (Estimate - TrueEffect) / (TrueEffect),
    parameter = factor(parameter, 
                       levels = c("pfda","pfna","pfba","lbpfoa","lbpfos","lbpfhxs","pfos")),
    RelativeCIwidth = ifelse(TrueEffect == 0, NA, CI_width / abs(TrueEffect))
  )


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


coefPFASstrong <- subset(coefPFAS, parameter == "pfda"| parameter =="pfna")
name_map <- c("pfda"="PFDA","pfna"="PFNA")
coefPFASstrong$parameter<-name_map[coefPFASstrong$parameter]
coefPFASstrong$parameter<-factor(coefPFASstrong$parameter , name_map)


### Plot absolute bias
#---------------------
p_joint_bias <- coefPFASstrong %>%
  filter(Group!="Semi-parametric")%>%
  ggplot( aes(x = parameter, y = RelativeBias, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.2, alpha = 0.1, show.legend = FALSE) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 16),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) +
  guides(fill = guide_legend(nrow = 4)) +  # Legend with 4 rows
  ggtitle("Relative bias of the individual effect estimate ") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") +
  labs(fill = 'Method:')


### Plot the CI width
#--------------------
p_joint_width <- coefPFASstrong %>%
  filter(Group!="Semi-parametric" & Method!="WQS (one) regression (rh=1 , bs=100)" &
           Method!="WQS (abst) regression (rh=1 , bs=100)" &
           Method!="WQS (one) regression (100 random subsets)")%>%
  ggplot( aes(x=parameter, y=RelativeCIwidth, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.2, alpha = 0.1, show.legend = FALSE) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 16),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) +
  guides(fill = guide_legend(nrow = 4)) +  
  labs(fill = 'Method:')+
  ggtitle("Relative width of 95% confidence/credible interval for the individual effect") +
  xlab("") + ylab("Interval width")


### Summarize empirical power using exact binomial confidence intervals
#----------------------------------------------------------------------
power_summary <- coefPFASstrong %>%
  group_by(Method, Size,Group,parameter) %>%
  summarise(
    Positives = sum(Power == 1, na.rm = TRUE),  # Count of significant results
    n = n(),                                    # Total number of simulations
    Power = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )


### PLot the empirical power
#---------------------------
p_joint_power <- power_summary %>%
  filter(Group!="Semi-parametric" & Group!="Weighted index")%>%
  filter(Method!= "WQS (one) regression (rh=1 ,bs=100") %>%
  ggplot( aes(x=parameter, y=Power, fill=Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
  ylim(0,1)+
  theme(
    plot.title = element_text(size = 16),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) +
  guides(fill = guide_legend(nrow = 4)) +  # Legend with 4 rows
  ggtitle("Empirical power to detect the individual effect") +
  ylab("Power") + xlab("")  +
  labs(fill = 'Method:')


### Summarize empirical CI coverage using exact binomial confidence intervals
#----------------------------------------------------------------------------
CI_summary <- coefPFASstrong %>%
  group_by(Method, Size, Group, parameter) %>%
  summarise(
    Positives = sum(CI_coverage == 1, na.rm = TRUE),  # Count of significant results
    n = n(),                                    # Total number of simulations
    CI_coverage = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )


### Plot the CI coverage
#-----------------------
p_joint_coverage <- CI_summary %>%
  filter(Group!="Semi-parametric" & Method!="WQS (one) regression (rh=1 , bs=100)" &
           Method!="WQS (abst) regression (rh=1 , bs=100)" &
           Method!="WQS (one) regression (100 random subsets)")%>%
  ggplot( aes(x=parameter, y=CI_coverage, fill=Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 16),
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) +
  guides(fill = guide_legend(nrow = 4)) +
  ylab("Coverage") + xlab("") +
  ggtitle("Empirical 95% confidence/credible interval coverage for the individual effect") +
  labs(fill = 'Method:')+
  geom_hline(yintercept = 0.95, col = "red") 



### Combine plots with a common legend
#-------------------------------------
final_plot <- ggarrange(
  p_joint_bias + theme(legend.position = "none"),
  p_joint_width + theme(legend.position = "none"),
  p_joint_coverage + theme(legend.position = "none"),
  p_joint_power + theme(legend.position = "none"),
  ncol = 2, nrow = 2,
  widths = c(0.5, 0.5), heights = c(0.6, 0.4),
  common.legend = TRUE,     
  legend = "bottom")


final_plot <- annotate_figure(
  final_plot,
  top = text_grob(
    "Simulation results for correlated predictors PFDA & PFNA (50% & 50% of total joint effect) – Individual effect estimation ", 
    face = "bold", size = 16
  )
)

final_plot



#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


name_map <- c("pfda"="PFDA (correlated + effect on outcome)","pfna"="PFNA (correlated + effect on outcome)",
              "pfba" = "PFBA (nuisance)","lbpfoa" = "PFOA (total) (nuisance)","lbpfos" = "PFOS (total) (nuisance)",
              "lbpfhxs" = "PFHXS (total) (nuisance)",
              "pfos" = "PFOS (nuisance)")
coefPFAS$parameter<-name_map[coefPFAS$parameter]
coefPFAS$parameter<-factor(coefPFAS$parameter , name_map)


final_plot <- coefPFAS %>%
  filter(Group!="OLS methods" & Group!="Weighted index"& Method!="Ridge regression (bs=200)" &
           Method!="Bayesian Lasso regression"& Method!="Bayesian horseshoe regression")%>%
  ggplot( aes(x = parameter, y = PIP, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 2) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 15),        # Remove x-axis labels
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2)) +  # Legend with 4 rows
  ylab("Posterior inclusion probability (PIP) or bootstrap inclusion probability (BIP)") +
  labs(fill = 'Method: ')


final_plot <- annotate_figure(
  final_plot,
  top = text_grob(
    "Inclusion probabilities for two correlated equal effects and 5 additional nuisance parameters", 
    face = "bold", size = 16
  )
)

final_plot

