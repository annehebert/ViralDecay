
library(dplyr)
library(ggplot2)

# load data files
data = read.csv("pave_data_formatted.csv") # data points

indFilePath="PAVE_biphasic_decay_final/IndividualParameters/estimatedIndividualParameters.txt" 
indParams = read.csv(indFilePath) #individual fit params

popFilePath="tabestimates_best_finalrun.csv"
popParams = read.csv(popFilePath) #population params

# define double exp
double.exp <- function(v0, a, b1, b2, t){
  y = v0*((1-a)*exp(-b1*t*7)+a*exp(-b2*t*7))
  return(y)
}

# prepare individual fits
format.ind.data <- function(estimate='mode'){
  # set estimate to either 'mode', 'mean', or 'SAEM'
  # extracts individual parameters
  # evaluates resulting double exponential
  
  # define time vector over which to evaluate double exponential
  t_vec <- seq(0,48,0.1)
  # create empty dataframe
  df <- NULL
  df <- data.frame(ids=character(), group=numeric(), 
                   vals=numeric(), weeks=numeric())
  
  # extract individual fit params
  for (ind in indParams$id){
    # set individual params according to estimate type
    tempParams = indParams[indParams$id==ind,]
    if (estimate=='mode'){
      x0_temp = 10^(tempParams$x0_mode)
      a_temp = 10^(tempParams$a_mode)
      b1_temp = tempParams$b1_mode
      b2_temp = tempParams$b2_mode
    } else if (estimate=='SAEM'){
      x0_temp = 10^(tempParams$x0_SAEM)
      a_temp = 10^(tempParams$a_SAEM)
      b1_temp = tempParams$b1_SAEM
      b2_temp = tempParams$b2_SAEM
    } else if (estimate=='mean'){
      x0_temp = 10^(tempParams$x0_mean)
      a_temp = 10^(tempParams$a_mean)
      b1_temp = tempParams$b1_mean
      b2_temp = tempParams$b2_mean
    }
    
    # evaluate double exp with individual params and create dataframe
    vals = double.exp(x0_temp, a_temp, b1_temp, b2_temp,t_vec)
    ids = rep(ind, length(vals))
    group = rep(tempParams$group, length(vals))
    weeks = t_vec
    indFits_temp = data.frame(ids, vals, group, weeks)
    df = rbind(df, indFits_temp)
  }
  return(df)
}

indFits <- format.ind.data(estimate ='mode')

# INDIVIDUAL FITS
# turn variables into factors as needed. sort data by group, and then by animal ID, making sure that the animal ID variable encodes this ordering on it's own so that when we plot it's in the right order 
indFits <- indFits %>% mutate(group = as.factor(group), ID = as.factor(ids)) %>%  
  arrange(group,ID) %>%
  mutate(group_ID = as.factor(paste(group,ID,sep="_")))

# DATA POINTS
# turn variables into factors as needed. sort data by group, and then by animal ID, making sure that the animal ID variable encodes this ordering on it's own so that when we plot it's in the right order 
data <- data %>% mutate(group = as.factor(group), ID = as.factor(ID),censor = as.factor(censor),RhmAbs = as.factor(RhmAbs), N803 = as.factor(N803)) %>%  
  arrange(group,ID) %>%
  mutate(group_ID = as.factor(paste(group,ID,sep="_")))


# extracting the order of the animals (by ID) when they are ordered by group then ID
ID_levels <- sapply(strsplit(levels(data$group_ID),"_"),'[',2)

# applying this ordering to the ID variable for data points and individual fits
data <- data %>%
  mutate(ID = factor(ID,levels=ID_levels))
indFits <- indFits %>%
  mutate(ID = factor(ID,levels=ID_levels))


ggplot() + 
  #uncensored data
  geom_point(data = data, aes(x = Weeks, y = level, color = group, shape = censor)) + 
  # individual fits
  geom_line(data = indFits, aes(x = weeks, y = vals)) +
  facet_wrap(
    ~ID,
    ncol = 7,
    scales = "fixed",
    axes = "all",
    axis.labels = "all_x"
  ) + 
  geom_hline(yintercept = 59, linetype = "dashed") +
  labs(x = 'Time after treatment initiation (weeks)', y = 'Plasma SIV RNA\n (copies/mL)') +
  theme_classic() +
  scale_y_log10(
    breaks = scales::trans_breaks("log10",function(x) 10^x),
    labels = scales::trans_format("log10",scales::math_format(10^.x)),
    limits = c(1e-4,1e8)
  ) + 
  scale_x_continuous(
    limits = c(0,30),
    minor_breaks = seq(0,30,5),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_shape_manual(values=c(19,21)) + 
  theme(plot.title = element_text(hjust = 0.5),strip.background = element_rect(color=NA))





# Population average plot

# extract parameters from best run
V0 <- 10^(popParams$estimate[popParams$param=='x0_pop'])
A <- 10^(popParams$estimate[popParams$param=='a_pop'])
b1 <- popParams$estimate[popParams$param=='b1_pop']
b2 <- popParams$estimate[popParams$param=='b2_pop']
beta_b1_N803 <- popParams$estimate[popParams$param=='beta_b1_N803_1']
beta_b1_RhmAbs <- popParams$estimate[popParams$param=='beta_b1_RhmAbs_1']
beta_b2_N803 <- popParams$estimate[popParams$param=='beta_b2_N803_1']
time_vec <- seq(0,30,0.05)

# evaluate exp trajectory for each group
group1 <- double.exp(V0, A, b1, b2, time_vec)  
group2 <- double.exp(V0, A, b1+beta_b1_RhmAbs, b2, time_vec)  
group3 <- double.exp(V0, A, b1+beta_b1_N803+beta_b1_RhmAbs, b2+beta_b2_N803, time_vec) 

# dataframe containing trajectory for each group
data_groups <- data.frame(time_vec, group1, group2, group3)

ggplot(data=data_groups)+
  geom_line(aes(x = time_vec,  y = group1), linewidth=1)+
  geom_line(aes(x = time_vec,  y = group2), linewidth=1)+
  geom_line(aes(x = time_vec,  y = group3), linewidth=1)+
  theme_classic()+
  geom_hline(yintercept=59, linetype = "dashed")+
  xlab('Time after treatment initiation (weeks)')+ylab('Plasma SIV RNA \n(copies/mL)')+
  scale_x_continuous(
    limits = c(0,30),
    minor_breaks = seq(0,30,5),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10",function(x) 10^x),
    labels = scales::trans_format("log10",scales::math_format(10^.x)),
    limits = c(1e-4,1e8)
  )


