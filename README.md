# ViralDecayFits

Contains Monolix code for modeling viral decay upon treatment initiation with a biphasic exponential decay. R script for performing repeated parameter estimation once statistical and covariate model have been determined, and for generating some plots.

## Files

- batch_runs.R : R code for running repeated parameter estimations while varying initial parameter values, and for generating some plots
- PAVE_biphasic_decay_final.mlxtran : Monolix file for initial analysis and model selection
- exp_decay_model.txt : file defining the biphasic decay model (log10 transformed)
