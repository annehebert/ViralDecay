
library(dplyr)
library(ggplot2)

data = read.csv("pave_data_formatted.csv")

indFilePath="PAVE_biphasic_decay_final/IndividualParameters/estimatedIndividualParameters.txt" 

indParams = read.csv(indFilePath)

# turn variables into factors as needed. sort data by group, and then by animal ID, making sure that the animal ID variable encodes this ordering on it's own so that when we plot it's in the right order 
data <- data %>% mutate(group = as.factor(group), ID = as.factor(ID),censor = as.factor(censor),RhmAbs = as.factor(RhmAbs), N803 = as.factor(N803)) %>%  
  arrange(group,ID) %>%
  mutate(group_ID = as.factor(paste(group,ID,sep="_")))

# extracting the order of the animals (by ID) when they are ordered by group then ID
ID_levels <- sapply(strsplit(levels(data$group_ID),"_"),'[',2)

# applying this ordering to the ID variable
data <- data %>%
  mutate(ID = factor(ID,levels=ID_levels))

ggplot() + 
  #uncensored data
  geom_point(data = data, aes(x = Weeks, y = level, color = group, shape = censor)) + 
  facet_wrap(
    ~ID,
    ncol = 7,
    scales = "free",
    axis = "all"
  ) + 
  labs(x = 'Time after treatment initiation (weeks)', y = 'Plasma SIV RNA\n (copies/mL)', title = data$ID[1]) +
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



