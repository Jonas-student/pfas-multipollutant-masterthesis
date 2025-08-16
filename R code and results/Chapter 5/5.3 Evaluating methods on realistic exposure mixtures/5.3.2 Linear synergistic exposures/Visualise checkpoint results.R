################################################################################
# Simulation study Master Thesis - Jonas Meijerink
# Linear synergestic exposures visualisations
# 27/06/2025
################################################################################


library(ggpubr)
library(hrbrthemes)
library(viridis)
library(ggplot2)


### Load previous results
#-------------------------
output <- do.call(rbind, checkpoint_results$results)


################################################################################

#############################  joint effect  ###################################

################################################################################


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


### Plot absolute bias
#---------------------
p_joint_bias <- ggplot(coefwqs, aes(x = Method, y = Estimate, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.2, alpha = 0.1, show.legend = FALSE) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
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
  ) +
  guides(fill = guide_legend(nrow = 4)) +  # Legend with 4 rows
  ggtitle("Absolute bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") +
  labs(fill = 'Method:')


### Plot the CI width
#--------------------
p_joint_width <- coefwqs %>%
  filter(Method != "Bayesian model averaging") %>%
  ggplot( aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.2, alpha = 0.1, show.legend = FALSE) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
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
  ) +
  guides(fill = guide_legend(nrow = 4)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Summarize empirical power using exact binomial confidence intervals
#----------------------------------------------------------------------
power_summary <- coefwqs %>%
  filter(Method != "Bayesian model averaging") %>%
  group_by(Method, Size,Group) %>%
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
  ggplot( aes(x=Method, y=Power, fill=Method)) +
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
    axis.text.x = element_blank(),        # Remove x-axis labels
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",           # Move legend to bottom
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11.5)
  ) +
  guides(fill = guide_legend(nrow = 4)) +  # Legend with 4 rows
  ggtitle("Empirical power") +
  ylab("Power") + xlab("")  +
  labs(fill = 'Method:')


#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#


### Summarize empirical CI coverage using exact binomial confidence intervals
#----------------------------------------------------------------------------
CI_summary <- coefwqs %>%
  filter(Method != "Bayesian model averaging") %>%
  group_by(Method, Size, Group) %>%
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
  ggplot( aes(x=Method, y=CI_coverage, fill=Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  theme_ipsum() +
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
  ) +
  guides(fill = guide_legend(nrow = 4)) +
  ylab("Coverage") + xlab("") +
  ggtitle("Empirical 95% confidence/credible interval coverage") +
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
  widths = c(0.54, 0.46), heights = c(0.6, 0.4),
  common.legend = TRUE,     
  legend = "bottom")
final_plot <- annotate_figure(
  final_plot,
  top = text_grob(
    "Simulation results linear synergistic exposures - Joint effect estimation", 
    face = "bold", size = 16
  )
)

final_plot
