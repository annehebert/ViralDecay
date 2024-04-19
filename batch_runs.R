# A simple script for running the parameter estimation for a simple biphasic 
# exponential decay model with random effects and two categorical covariates

# I wrote this to use once the model selection had been finalized using the
# Stepwise Covariate Modeling algorithm in Monolix's "model building" section

# credit to monolix website, I copied heavily from their examples

# load libraries and initialize monolix API
library(ggplot2)
library(gridExtra)
library(lixoftConnectors)
library(Rsmlx)
library(dplyr)
initializeLixoftConnectors(software="monolix")

# load project
project <- paste0("PAVE_allgroups_run21.mlxtran")
loadProject(projectFile = project)
# double check model is defined correctly
getIndividualParameterModel()
# if needed, set the covariate model
#setCovariateModel(b1=c(N803=TRUE, RhmAbs=TRUE), b2=c(N803=TRUE))


# ----------------------------------------------------------------------------
# Run parameter estimations with resampled initial estimates,
# store the parameter estimation results in tabestimates
# store the p values (Wald test) in pvals
# ----------------------------------------------------------------------------

popparams <- getPopulationParameterInformation()
setPopulationParameterEstimationSettings(nbBurningIterations=5,
                                         exploratoryInterval=200,
                                         smoothingInterval=100,
                                         smoothingAlpha=0.7)
#setInitialEstimatesToLastEstimates()
estimates <- NULL;
tabestimates <- NULL; tabiters <- NULL; pvals <- NULL;

for(i in 1:100){
  # sample new initial estimates
  popini <- sapply(1:nrow(popparams), 
                   function(j){
                     if(popparams$initialValue[j]>=0){
                       runif(n=1, 
                             min=max(0,popparams$initialValue[j]/2), 
                             max=popparams$initialValue[j]*2
                             )
                     }else{-runif(n=1, 
                                  min=-popparams$initialValue[j]/2, 
                                  max=-popparams$initialValue[j]*2
                                  )
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
  
  # store the covariate parameter p values
  pvals_temp <- as.data.frame(getTests()$wald$p.value)
  colnames(pvals_temp) <- 'pvalue'
  pvals_temp$parameters <- getTests()$wald$parameter
  pvals_temp$run <- i
  pvals <- rbind(pvals, pvals_temp)
  
  # store the iterations
  iters <- getChartsData("plotSaem")
  iters$run <- i
  tabiters <- rbind(tabiters, iters)
}



# plot SAEM iterations
plotList <- list()
i <- 1
for (param in popparams$name) {
  if (popparams[popparams$name == param, ]$method == "FIXED") next
  changePhase <- tabiters$iteration[which(diff(tabiters$phase) == 1) + 1]
  plotList[[i]] <- ggplot(tabiters, aes_string(x = "iteration", y = param)) +
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
  estimates <- tabestimates[tabestimates$param == param, ]
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

#getTests()
#write.csv(tabiters, 'tabiters_finalrun.csv')  


# find best run
tabestimates_best <- tabestimates %>%
  arrange(ofv)
tabestimates_best <- tabestimates_best[1:11,]
best_run <- tabestimates_best$run[1]
pvals_best <- pvals %>% filter(run==14)

# save parameters from best run
x0 <- tabestimates_best$estimate[tabestimates_best$param=='x0_pop']
a <- tabestimates_best$estimate[tabestimates_best$param=='a_pop']
b1 <- tabestimates_best$estimate[tabestimates_best$param=='b1_pop']
b2 <- tabestimates_best$estimate[tabestimates_best$param=='b2_pop']
beta_b1_N803 <- tabestimates_best$estimate[tabestimates_best$param=='beta_b1_N803_1']
beta_b1_RhmAbs <- tabestimates_best$estimate[tabestimates_best$param=='beta_b1_RhmAbs_1']
beta_b2_N803 <- tabestimates_best$estimate[tabestimates_best$param=='beta_b2_N803_1']

# save RSEs from best run
x0_se <- x0/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='x0_pop'])
a_se <- a/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='a_pop'])
b1_se <- b1/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='b1_pop'])
b2_se <- b2/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='b2_pop'])
beta_b1_N803_se <-beta_b1_N803/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='beta_b1_N803_1'])
beta_b1_RhmAbs_se <-beta_b1_RhmAbs/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='beta_b1_RhmAbs_1'])
beta_b2_N803_se <- beta_b2_N803/100*as.numeric(tabestimates_best$rses[tabestimates_best$param=='beta_b2_N803_1'])

# define function for calculating 95% CIs from the param estimates and standard errors
calculate.confint <- function(param, se){
  upper = param+1.96*se
  lower = param-1.96*se
  return(c(lower, upper))
}


# standard errors for additive treatment effects
b1_betaRhmAbs_se <- sqrt(b1_se^2+beta_b1_RhmAbs_se^2)
b1_betaRhmAbsN803_se <- sqrt(b1_betaRhmAbs_se^2+beta_b1_N803_se^2)
b2_betaN803_se <- sqrt(b2_se^2+beta_b2_N803_se^2)

calculate.confint(x0, x0_se)
calculate.confint(a, a_se)
calculate.confint(b1, b1_se)
calculate.confint(b2, b2_se)
calculate.confint(beta_b1_N803, beta_b1_N803_se)
calculate.confint(beta_b1_RhmAbs, beta_b1_RhmAbs_se)
calculate.confint(beta_b2_N803, beta_b2_N803_se)

calculate.confint(param = b1 + beta_b1_RhmAbs, 
                  se = b1_betaRhmAbs_se)
calculate.confint(param = b1+beta_b1_N803+beta_b1_RhmAbs, 
                  se = b1_betaRhmAbsN803_se)
calculate.confint(param = b2+beta_b2_N803,
                  se = b2_betaN803_se)




## Generate plots ------

# plot all raw data
plotObservedData()

# plot each individual fit, by group
plotIndividualFits(settings = list(cens=FALSE,
                                   grid=FALSE,
                                   ncol=7,
                                   xlim=c(0,30),ylim=c(-5,10)),
                   stratify = list(colorGroup = list(list(name = 'group')))
)+geom_hline(yintercept=log10(60))


# plot mean trajectory for each group
double.exp <- function(v0, a, b1, b2, t){
  y = v0 + log10((1-10^a)*exp(-b1*t*7) + 10^a*exp(-b2*t*7))
  return(y)
}


time_vec <- seq(0,30,0.05)
group1 <- double.exp(x0, a, b1, b2, time_vec)  
group2 <- double.exp(x0, a, b1+beta_b1_RhmAbs, b2, time_vec)  
group3 <- double.exp(x0, a, b1+beta_b1_N803+beta_b1_RhmAbs, b2+beta_b2_N803, time_vec)  

# dataframe containing trajectory for each group
data_groups <- data.frame(time_vec, group1, group2, group3)
  
ggplot(data=data_groups)+
  geom_line(aes(x = time_vec,  y = group1), linewidth=1)+
  geom_line(aes(x = time_vec,  y = group2), linewidth=1)+
  geom_line(aes(x = time_vec,  y = group3), linewidth=1)+
  theme_classic()+
  geom_hline(yintercept=log10(59))+
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3,4,5,6,7),
                    labels=c('','-2','','0','','2','','4','','6',''))
  

# plot estimated group mean and. individual parameters by group 
x0_estimates <- getEstimatedIndividualParameters()$saem$x0
a_estimates <- getEstimatedIndividualParameters()$saem$a
b1_estimates <- getEstimatedIndividualParameters()$saem$b1
ids <- getEstimatedIndividualParameters()$saem$id
group <- rep(1, length(ids))
group[8:14] <- 2
group[15:19] <- 3

param_estimates_df <- data.frame(group, ids, x0_estimates, a_estimates, b1_estimates)

param_means <- param_estimates_df %>%
  group_by(group) %>%
  summarize(x0_means = mean(x0_estimates),
            a_means = mean(a_estimates),
            b1_means = mean(b1_estimates))

p1 <- ggplot()+
  geom_jitter(data=param_estimates_df,
             aes(x = as.factor(group), y = x0_estimates, 
                 color=as.factor(group)),
             width = 0.25)+
  geom_hline(yintercept = param_means$x0_means[1])+
  geom_hline(yintercept = param_means$x0_means[2])+
  geom_hline(yintercept = param_means$x0_means[3])+
  theme_bw()+
  theme(legend.position = 'none')+
  ylim(4,8)

p2 <- ggplot()+
  geom_jitter(data=param_estimates_df,
              aes(x = as.factor(group), y = a_estimates, 
                  color=as.factor(group)),
              width = 0.25)+
  geom_hline(yintercept = param_means$a_means[1])+
  geom_hline(yintercept = param_means$a_means[2])+
  geom_hline(yintercept = param_means$a_means[3])+
  theme_bw()+
  theme(legend.position = 'none')+
  ylim(-6,-2)

p3 <- ggplot()+
  geom_jitter(data=param_estimates_df,
              aes(x = as.factor(group), y = b1_estimates, 
                  color=as.factor(group)),
              width = 0.25)+
  geom_hline(yintercept = param_means$b1_means[1])+
  geom_hline(yintercept = param_means$b1_means[2])+
  geom_hline(yintercept = param_means$b1_means[3])+
  theme_bw()+
  theme(legend.position = 'none')+
  ylim(0.5,1.5)

grid.arrange(p1,p2,p3, ncol=1)
