# Section 5.3 — Evaluating methods on realistic exposure mixtures

This folder contains all the code and resources used for the simulation study in **Section 5.3** of the thesis.

## Structure

### 1. R scripts ending on `_store_results.R`
The R scripts provide the code to analyse a dataset with a specific method. 
The input is a simulated dataset and the ouput is data frame containting all the estimates and metadata.

### 2. Subfolders for each simulation study
This folder contains:
- The `... simulation procedure.R` used to generate the datasets  
- The `checkpoint_results.rds` contains all estimated parameter and metadata 
- The `Visualise checkpoint results.R` 

- **Figures** describing the results
