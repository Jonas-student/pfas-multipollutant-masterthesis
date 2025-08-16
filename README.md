# Evaluating multi-pollutant methods for PFAS mixture effects on immunometabolic health – Master's Thesis

This repository contains the R code and simulation results used in my master’s thesis:  
**"Evaluating multi-pollutant methods for PFAS mixture effects on immunometabolic health"**  
This master's thesis was the result of a collaboration between Hasselt University and VITO (Flemish Institute for Technological Research).



## Repository Structure

The repository is organised by thesis chapters outlined below. The data used in this project originates from VITO (Flemish Institute for Technological Research) and is part of the ``Teenager HBM Study - 3M site''. Due to privacy regulations and data protection agreements, this dataset is not publicly available and is not included in this repository. However, the code, workflow and analysis pipeline provided here can help others understand how different multipollutant methods were applied in this thesis. 


## Contents Overview

### Chapter 1. Introduction
- Background  
- Societal relevance and stakeholder awareness  
- Research question  
- Ethical considerations  
- Structure of the thesis  

### Chapter 2. Data
- Human biomonitoring studies
- [Case study](/R%20code%20and%20results/Chapter%202/)

### Chapter 3. Shrinkage methods
- Frequentist shrinkage methods: Ridge, LASSO, Elastic Net  
- Bayesian shrinkage methods: Bayesian Lasso, spike-and-slab, horseshoe priors  

### Chapter 4. Multi-pollutant Methods
- Weighted Quantile Sum (WQS) regression and its variants  
- G-computation and causal inference framework  
- Bayesian adaptive sampling for variable selection and model averaging  
- Bayesian kernel machine regression (BKMR)  

### Chapter 5. Simulation Results
- [Continuous vs quantised exposures in WQS](/R%20code%20and%20results/Chapter%205/5.1%20Continuous%20vs%20quantised%20exposures%20in%20WQS/)  
- [Joint effect estimation using OLS](/R%20code%20and%20results/Chapter%205/5.2%20Joint%20effect%20estimation%20using%20OLS/)
- Evaluating methods on realistic exposure mixtures: 
    - [linear additive exposures](/R%20code%20and%20results/Chapter%205/5.3%20Evaluating%20methods%20on%20realistic%20exposure%20mixtures/5.3.1%20Linear%20additive%20exposures.R/)
    - [linear synergistic exposures](/R%20code%20and%20results/Chapter%205/5.3%20Evaluating%20methods%20on%20realistic%20exposure%20mixtures/5.3.2%20Linear%20synergistic%20exposures/)
    - [evaluating grouping behaviour](/R%20code%20and%20results/Chapter%205/5.3%20Evaluating%20methods%20on%20realistic%20exposure%20mixtures/5.3.3%20Evaluating%20grouping%20behaviour/)

### Chapter 6. Case Study: PFAS exposure and immunometabolic health 
- [Leukocyte count](/R%20code%20and%20results/Chapter%206/Leukocyte%20count/)
- [CD4+/CD8+ T-cell ratio](/R%20code%20and%20results/Chapter%206/CD4+CD8+%20T-cel%20ratio/)

### Chapter 7. Discussion  
- Interpretation of results, methodological assumptions, limitations and directions for future research

### Chapter 8. Conclusion
- Summary of key findings

### Appendices
- Appendix A: Software details & AI tools
- Appendix B: Causal assumptions
- [Appendix C: Additional simulation results and insights](/R%20code%20and%20results/Appendix%20C/)
- [Appendix D: Details on simulation procedure](/R%20code%20and%20results/Appendix%20D/)
- [Appendix E: Case study: Descriptive statistics, results and implementation](/R%20code%20and%20results/Appendix%20E/)



