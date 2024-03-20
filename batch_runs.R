# A simple script for running the parameter estimation for a simple biphasic 
# exponential decay model with random effects and two categorical covariates

# I wrote this to use once the model selection had been finalized using the
# Select Covariate Modeling algorithm in Monolix's "model building" section

# credit to monolix website, I copied heavily from one of their examples

# load libraries and initialize monolix API
library(ggplot2)
library(gridExtra)
library(lixoftConnectors)
library(dplyr)
initializeLixoftConnectors(software="monolix")

# load project
project <- paste0("PAVE_allgroups_run21.mlxtran")
loadProject(projectFile = project)
# double check model is defined correctly
getIndividualParameterModel()
setCovariateModel(b1=c(N803=TRUE, RhmAbs=TRUE))



# ----------------------------------------------------------------------------
# Run parameter estimations with different initial estimates,
# store the results in tabestimates
# ----------------------------------------------------------------------------
popparams <- getPopulationParameterInformation()
setPopulationParameterEstimationSettings(nbBurningIterations=5,
                                         exploratoryInterval=100,
                                         smoothingInterval=50,
                                         smoothingAlpha=0.7)
tabestimates <- NULL; tabiters <- NULL
for(i in 1:100){
  # sample new initial estimates
  popini <- sapply(1:nrow(popparams), 
                   function(j){
                     if(popparams$initialValue[j]>=0){runif(n=1, 
                                     min=max(0,popparams$initialValue[j]/2), 
                                     max=popparams$initialValue[j]*2)
                       }else{-runif(n=1, 
                               min=-popparams$initialValue[j]/2, 
                               max=-popparams$initialValue[j]*2)
                       }
                     })
  
  # set sampled values as new initial estimates
  newpopparams <- popparams
  newpopparams$initialValue <- popini
  setPopulationParameterInformation(newpopparams)
  
  # run the estimation
  runScenario()
  
  # store the estimates and s.e. in a table
  estimates <- as.data.frame(getEstimatedPopulationParameters())
  names(estimates) <- "estimate"
  rses <- getEstimatedStandardErrors()$stochasticApproximation$rse
  names(rses) <- getEstimatedStandardErrors()$stochasticApproximation$parameter
  rses <- as.data.frame(rses)
  estimates <- merge(estimates, rses, by = "row.names")
  estimates$run <- i
  estimates$ofv <- getEstimatedLogLikelihood()$importanceSampling[1]
  estimates$bic <- getEstimatedLogLikelihood()$importanceSampling[3]
  names(estimates)[names(estimates) == "Row.names"] <- "param"
  tabestimates <- rbind(tabestimates, estimates)
  
  # store the iterations
  iters <- getChartsData("plotSaem")
  iters$run <- i
  tabiters <- rbind(tabiters, iters)
}

tabestimates_filtered_with_bic <- tabestimates %>%
  filter(ofv < 350)
selected_runs <- unique(tabestimates_filtered$run)
tabiters_filtered <- tabiters %>%
  filter(run %in% selected_runs)

# plot SAEM iterations
plotList <- list()
i <- 1
for (param in popparams$name) {
  if (popparams[popparams$name == param, ]$method == "FIXED") next
  changePhase <- tabiters_filtered$iteration[which(diff(tabiters_filtered$phase) == 1) + 1]
  plotList[[i]] <- ggplot(tabiters_filtered, aes_string(x = "iteration", y = param)) +
    geom_line(aes(group = run, color = factor(run))) +
    theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
    geom_vline(xintercept = changePhase, color = 1:length(changePhase)) +
    labs(title = param, x = NULL, y = NULL)
  i <- i + 1
}
grid.arrange(grobs = plotList, ncol = 3)

# plot population parameters
plotList <- list()
i <- 1
for (param in popparams$name) {
  if (popparams[popparams$name == param, ]$method == "FIXED") next
  estimates <- tabestimates_filtered[tabestimates_filtered$param == param, ]
  plotList[[i]] <- ggplot(estimates, aes(x = run, y = estimate)) +
    geom_point(aes(color = factor(run))) +
    geom_errorbar(aes(ymax = estimate * (1+as.numeric(rses)/100), 
                      ymin = estimate * (1-as.numeric(rses)/100), 
                      color = factor(run))) +
    theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
    labs(title = param, x = NULL, y = NULL)
  i <- i + 1
}
grid.arrange(grobs = plotList, ncol = 3)

param_means <- tabestimates_filtered %>%
  group_by(param) %>%
  summarize(mean = mean(estimate))


#write.csv(tabestimates, 'param_estimates.csv')  

tabestimates_best <- tabestimates %>%
  arrange(ofv)
tabestimates_best <- tabestimates_best[1:11,]

