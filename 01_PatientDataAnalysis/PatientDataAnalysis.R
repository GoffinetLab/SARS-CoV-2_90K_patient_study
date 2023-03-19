## load required libraries ##

library(ggplot2)
library(tidyverse)
library(ggfortify)
library(gridExtra)
library(nlme)
library(MuMIn)
library(predictmeans)
library(patchwork)
library(DT)
library(emmeans)
library(ggpubr)
library(wesanderson)

#### Read in data ####

alldat <- read.csv("path/to/90K_data_all.csv",
                   sep=";")

#### Figure 1 ####

## 1A

# rework data to include only samples with matched healthy controls

# create separate entries for healthy and covid19 patients
sick <- alldat %>%
  #mutate(serum_90k = as.numeric(serum_90k)) %>%
  subset(!is.na(serum_90k)) %>%
  select(pat_id, dpso, serum_90k, age, sex, who) %>%
  mutate(condition = "covid19") %>%
  mutate(match = pat_id) %>%
  mutate(dpso = ifelse(dpso == "asymptomatisch", NA, dpso)) %>%
  mutate(dpso = as.numeric(dpso))

hc <- alldat %>%
  subset(!is.na(healthy_serum_90k)) %>%
  select(pat_id, dpso, healthy_serum_90k, age, sex, who) %>%
  rename(serum_90k = healthy_serum_90k) %>%
  mutate(condition = "healthy") %>%
  mutate(match = pat_id) %>%
  mutate(pat_id = paste0(pat_id, "_hc")) %>%
  mutate(dpso = 0)
# join
hc_match <- dplyr::full_join(sick, hc) %>%
  mutate(condition = fct_relevel(condition, c("healthy", "covid19")))


## do comparison using a linear mixed model

lme_hcmatch1 <- nlme::lme(serum_90k ~ condition, random = ~1|pat_id, data= hc_match)
summary(lme_hcmatch1)


## perform comparison of 90k levels (after log normalisation) between individual WHO levels or for grouping into moderate and severe diseas (3-4 vs 5-7)

## build model with who levels
# with each level individually
lme_who1 <- nlme::lme(serum_90k ~ who, random = ~1|pat_id, data= sick_who)
# moderate or severe disease
lme_who2 <- nlme::lme(serum_90k ~ who_level, random = ~1|pat_id, data= sick_who)

summary(lme_who1)
summary(lme_who2)

## compare using emmeans to estimate model means and perform pairwise contrasts
lme_who1.emm <- emmeans(lme_who1, "who")
lme_who2.emm <- emmeans(lme_who2, "who_level")


### 1C ###

## perform linear mixed effects model analysis to investigate whether 90k serum levels are different for different WHO severity levels of COVID-19 and if this changes over time

## log transform serum 90k level
alldat$log_serum_level <- log(alldat$serum_90k, base=10)
## assign WHO severity level
alldat$who_level <- ifelse(alldat$who %in% c(3:4), "3_4","5_7")
## REMOVE DROPOUT PATIENTS!
alldat <- alldat[alldat$pat_id != 30 & alldat$pat_id != 33,]

## build linear mixed effects model

model_lme_ad <- nlme::lme(log_serum_level ~ dpso + who_level, random = ~1|pat_id, data=alldat)
summary(model_lme_ad)

# graph result 

# plot predicted values based on model
newdat <- expand.grid(who_level=unique(alldat$who_level),
                      dpso=c(min(alldat$dpso),
                                        max(alldat$dpso)))
# make graph
ggplot(alldat, aes(x=days_post_onset, y=log_serum_level, colour=who_level)) +
  geom_point(size=3) +
  geom_line(data=newdat, aes(y=predict(model_lme_ad, level=0, newdata=newdat), size="Population")) +
  geom_line(aes(group=pat_id), linetype="dashed")+
  theme_classic() 


#### Figure 2 ####

### 2A ###

# subset alldat to get pbmc data

pbmc <- alldat %>%
  mutate(protein_pbmc_90k = as.numeric(protein_pbmc_90k)) %>%
  subset(!is.na(protein_pbmc_90k)) %>%
  mutate(condition = ifelse(grepl("hc", pat_id), "healthy", "covid19")) %>%
  mutate(condition = fct_relevel(condition, c("healthy", "covid19")))

## build linear mixed effects model comparing pbmc 90k levels COVID19 to healthy control

# log normalise pbmc 90k pbmc protein values
pbmc$log10_prot_pbmc <- log(pbmc$protein_pbmc_90k, base=10)
lme_pbmcvshc1 <- nlme::lme(log10_prot_pbmc ~ condition, random = ~1|pat_id, data= pbmc)
 summary(lme_pbmcvshc1)

 
## build linear mixed effects model comparing pbmc 90k levels COVID19 to healthy control

# assign WHO severity classification
pbmc_who <- pbmc %>%
   subset(condition == "covid19") %>%
   mutate(who_level = ifelse(who %in% c(3,4), "3_4", "5_7"))
# build model 
lme_pbmcwho1 <- nlme::lme(protein_pbmc_90k ~ who_level, random = ~1|pat_id, data= pbmc_who)
summary(lme_pbmcwho1)

#### Figure 3 ####

### 3A ###

## compare dCt values from PBMC qPCR data between COVID19 and healthy controls

#subset to get qpcr data only
qpcr <- alldat %>%
  mutate(qpcr_dct_90k = as.numeric(qpcr_dct_90k)) %>%
  subset(!is.na(qpcr_dct_90k)) %>%
  mutate(condition = ifelse(grepl("hc", pat_id), "healthy", "covid19")) %>%
  mutate(condition = fct_relevel(condition, c("healthy", "covid19")))
# build model
lme_pbmcrna1 <- nlme::lme(qpcr_dct_90k ~ condition, random = ~1|pat_id, data= qpcr)
summary(lme_pbmcrna1)


## compare dCt values from PBMC qPCR data between WHO severity levels

# assign WHO level
qpcr_who <- qpcr %>%
  subset(condition == "covid19") %>%
  mutate(who_level = ifelse(who %in% c(3,4), "3_4", "5_7"))
# build model
lme_pbmcrnawho1 <- nlme::lme(qpcr_dct_90k ~ who_level, random = ~1|pat_id, data= qpcr_who)
 summary(lme_pbmcrnawho1)

#### Sup. Figure 1 ####

## compare 90k serum concentrations in Dex treated vs non-Dex treated patients
 
# subset data to get data on dex treatment only
dex <- alldat %>%
   mutate(serum_90k = as.numeric(serum_90k)) %>%
   subset(!is.na(serum_90k) & !is.na(dex))
# build model
lme_dex1 <- nlme::lme(serum_90k ~ dex, random = ~1|pat_id, data= dex)
 summary(lme_dex1)
 