# BME_main
This repo contains codes for simulation study results for the manuscript "Bayesian Model Evaluation using $p$-values".

These codes were run on a computing cluster using slurm array jobs.

The codes are organized as follows:
1. files start with "nlr" are for the normal linear regression example, "gglm" for the gamma GLM example, and "sm" for the competing risk survival model example
2. For nlr runs:
  + run number 1 is used for level testing
  + run number 2-5 are for four alternative models for power testing
3. For gglm runs:
  + run number 1 is used for level testing for the conditional model
  + run number 2 is used for level testing for the joint model
  + run number 3-5 are for three alternative models for power testing
4. For sm runs:
  + run number 1 is used for level testing
  + run number 2-3 are for two alternative models for power testing
5. the `functions.R` file contains all the helper functions for simulation runs, aggregating results, and for plotting.
6. the `summary.R` file contains codes to aggregate run results as well as plotting.

Gamma GLM and competing risk survival model runs require `stan` files. They are both included here. In the R codes, these stan files were pre-compiled on the server and directly read in to save runtime. Users can obtain those pre-compiled files using the corresponding `stan` files. 
