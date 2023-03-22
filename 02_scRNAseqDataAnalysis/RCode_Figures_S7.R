library(lme4)
library(tidyverse)
library(broom.mixed)


#### Read in data ####

PBMC_scRNAseq_study_Berlin <- read.csv("path/to/PBMC_scRNAseq_study_Berlin.csv",sep=";")
PBMC_scRNAseq_study_Bonn <- read.csv("path/to/PBMC_scRNAseq_study_Bonn",sep=";")


### Mixed Linear Regression Model  ####

#### Figure Suppl. 7A ####


#### 1. Model with 2 groups (cases vs controls)

Classical_Monocytes_2g <- lmer(log(Classical_Monocytes+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLAmDRp_CD83p_Monocytes_2g <- lmer(log(HLAmDRp_CD83p_Monocytes+ 0.01) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

CD163p_Monocytes_2g <- lmer(log(CD163p_Monocytes+0.01) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLAmDRm_S100Ap_monocytes_2g <- lmer(log(HLAmDRm_S100Ap_monocytes+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

Non_classical_Monocytes_2g <- lmer(log(Non_classical_Monocytes+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)


dat_plot_2g <- data.frame(bind_rows(Classical_Monocytes_2g, HLAmDRp_CD83p_Monocytes_2g, CD163p_Monocytes_2g, HLAmDRm_S100Ap_monocytes_2g, Non_classical_Monocytes_2g)) %>% 
  filter(effect != "ran_pars", term != "(Intercept)") %>% 
  mutate(outcome = rep(c("log(Classical Monocytes)", "log(HLAmDRp CD83p Monocytes)", "log(CD163p Monocytes)", "log(HLADRm S100Ap Monocytes)", "log(Non-classical Monocytes)"), each=2))

dat_plot_2g_caco <- dat_plot_2g %>% 
  filter(term == "case_controlcase")
dat_plot_2g_caco$pos <- 1:nrow(dat_plot_2g_caco)
dat_plot <- dat_plot_2g_caco %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_2g_caco$text <- res1

ggplot(data=dat_plot_2g_caco, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-1.2, 10))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("difference between cases and control (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_2g_caco$outcome)), breaks=dat_plot_2g_caco$outcome)+
  geom_text(aes(x = pos, y = rep(8, each=nrow(dat_plot_2g_caco))), size=4.5,label = rev(res1))


dat_plot_2g_days <- dat_plot_2g %>% 
  filter(term == "days_post_symptoms")

dat_plot_2g_days$pos <- 1:nrow(dat_plot_2g_days)
dat_plot <- dat_plot_2g_days %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_2g_days$text <- res1

ggplot(data=dat_plot_2g_days, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-.5, 0.75))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("effect of days post symptom-onset (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_2g_days$outcome)), breaks=dat_plot_2g_days$outcome)+
  geom_text(aes(x = pos, y = rep(0.5, each=nrow(dat_plot_2g_days))), size=4.5,label = rev(res1))


#### 2. Model with 3 groups (control vs mild vs severe cases)

Classical_Monocytes_3g <- lmer(log(Classical_Monocytes+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLAmDRp_CD83p_Monocytes_3g <- lmer(log(HLAmDRp_CD83p_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

CD163p_Monocytes_3g <- lmer(log(CD163p_Monocytes+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLAmDRm_S100Ap_monocytes_3g <- lmer(log(HLAmDRm_S100Ap_monocytes+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value) 

Non_classical_Monocytes_3g <- lmer(log(Non_classical_Monocytes+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

dat_plot_3g <- data.frame(bind_rows(Classical_Monocytes_3g, HLAmDRp_CD83p_Monocytes_3g, CD163p_Monocytes_3g, HLAmDRm_S100Ap_monocytes_3g, Non_classical_Monocytes_3g)) %>% 
  filter(effect != "ran_pars", term != "(Intercept)") %>% 
  mutate(outcome = rep(c("log(Classical Monocytes)", "log(HLAmDRp CD83p Monocytes)", "log(CD163p Monocytes)", "log(HLA-DRm S100Ap Monocytes)", "log(Non-classical Monocytes)"), each=3))


dat_plot_3g_caco <- dat_plot_3g %>% 
  filter(term != "days_post_symptoms") %>% 
  mutate(group = rep(c("mild vs control", "severe vs control"), times=5))
dat_plot_3g_caco$pos <- c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2)
dat_plot <- dat_plot_3g_caco %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_3g_caco$text <- res1

ggplot(data=dat_plot_3g_caco, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high), col=rev(group))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-4, 14))+
  geom_hline(yintercept=0, lty=2) +
  scale_color_manual(values=cols[2:3])+
  xlab("") +
  ylab("difference between cases and control (95% CI)") +
  theme_bw(base_size = 20)+ 
  theme(legend.position=c(0.5,0.94), legend.title=element_blank(), legend.direction = "horizontal")+
  scale_x_discrete(limits = rev(c("", unique(dat_plot_3g_caco$outcome))), breaks=unique(dat_plot_3g_caco$outcome))+
  geom_text(aes(x = pos, y = rep(12, each=nrow(dat_plot_3g_caco))), size=4.5,label = rev(res1))


dat_plot_3g_days <- dat_plot_3g %>% 
  filter(term == "days_post_symptoms")

dat_plot_3g_days$pos <- 1:nrow(dat_plot_3g_days)
dat_plot <- dat_plot_3g_days %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_3g_days$text <- res1

ggplot(data=dat_plot_3g_days, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-.5, 0.75))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("effect of days post symptom-onset (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_2g_days$outcome)), breaks=dat_plot_2g_days$outcome)+
  geom_text(aes(x = pos, y = rep(0.5, each=nrow(dat_plot_2g_days))), size=4.5,label = rev(res1))




#### Figure Suppl. 7B ####


#### 1. Model with 2 groups (cases vs controls)

Classical_Monocytes_2g <- lmer(log(Classical_Monocytes+ 0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRhi_CD83hi_Monocytes_2g <- lmer(log(HLA_DRhi_CD83hi_Monocytes+ 0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRlo_CD163hi_Monocytes_2g <- lmer(log(HLA_DRlo_CD163hi_Monocytes+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRlo_S100Ahi_Monocytes_2g <- lmer(log(HLA_DRlo_S100Ahi_Monocytes+ 0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

Non_classical_Monocytes_2g <- lmer(log(Non_classical_Monocytes+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

mDCs_2g <- lmer(log(mDCs+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

pDCs_2g <- lmer(log(pDCs+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

NK_cells_2g <- lmer(log(NK_cells+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

Plasmablasts_2g <- lmer(log(Plasmablasts+0.001) ~ case_control + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)


dat_plot_2g <- data.frame(bind_rows(Classical_Monocytes_2g, HLA_DRhi_CD83hi_Monocytes_2g, HLA_DRlo_CD163hi_Monocytes_2g, HLA_DRlo_S100Ahi_Monocytes_2g, Non_classical_Monocytes_2g,
                                    mDCs_2g, pDCs_2g, NK_cells_2g, Plasmablasts_2g)) %>% 
  filter(effect != "ran_pars", term != "(Intercept)") %>% 
  mutate(outcome = rep(c("log(Classical Monocytes)", "log(HLA-DRhi CD83hi Monocytes)",  "log(HLA-DRlo CD163hi Monocytes)", "log(HLA-DRlo S100Ahi Monocytes)", "log(Non-classical Monocytes)", "log(mDCs)", "log(pDCs)", "log(NK cells)", "log(Plasmablasts)"), each=2))

dat_plot_2g_caco <- dat_plot_2g %>% 
  filter(term == "case_controlcase")
dat_plot_2g_caco$pos <- 1:nrow(dat_plot_2g_caco)
dat_plot <- dat_plot_2g_caco %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_2g_caco$text <- res1

ggplot(data=dat_plot_2g_caco, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-3, 15))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("difference between cases and control (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_2g_caco$outcome)), breaks=dat_plot_2g_caco$outcome)+
  geom_text(aes(x = pos, y = rep(12, each=nrow(dat_plot_2g_caco))), size=4.5,label = rev(res1))

dat_plot_2g_days <- dat_plot_2g %>% 
  filter(term == "days_post_symptoms")

dat_plot_2g_days$pos <- 1:nrow(dat_plot_2g_days)
dat_plot <- dat_plot_2g_days %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_2g_days$text <- res1

ggplot(data=dat_plot_2g_days, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-.5, 0.75))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("effect of days post symptom-onset (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_2g_days$outcome)), breaks=dat_plot_2g_days$outcome)+
  geom_text(aes(x = pos, y = rep(0.5, each=nrow(dat_plot_2g_days))), size=4.5,label = rev(res1))


#### 2. Model with 3 groups (control vs mild vs severe cases)

Classical_Monocytes_3g <-lmer(log(Classical_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRhi_CD83hi_Monocytes_3g <-lmer(log(HLA_DRhi_CD83hi_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRlo_CD163hi_Monocytes_3g <- lmer(log(HLA_DRlo_CD163hi_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

HLA_DRlo_S100Ahi_Monocytes_3g <- lmer(log(HLA_DRlo_S100Ahi_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

Non_classical_Monocytes_3g <- lmer(log(Non_classical_Monocytes+ 0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)


mDCs_3g <- lmer(log(mDCs+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

pDCs_3g <- lmer(log(pDCs+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

NK_cells_3g <- lmer(log(NK_cells+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

Plasmablasts_3g <- lmer(log(Plasmablasts+0.001) ~ disease_status + days_post_symptoms + (1|Patient_ID), data=PBMC_scRNAseq_study_Bonn) %>% 
  broom.mixed::tidy(conf.int=TRUE) %>% 
  dplyr::select(term, effect, group, estimate, conf.low, conf.high, p.value)

dat_plot_3g <- data.frame(bind_rows(Classical_Monocytes_3g, HLA_DRhi_CD83hi_Monocytes_3g, HLA_DRlo_CD163hi_Monocytes_3g, HLA_DRlo_S100Ahi_Monocytes_3g, Non_classical_Monocytes_3g,
                                    mDCs_3g, pDCs_3g, NK_cells_3g, Plasmablasts_3g)) %>% 
  filter(effect != "ran_pars", term != "(Intercept)") %>% 
  mutate(outcome = rep(c("log(Classical Monocytes)", "log(HLA-DRhi CD83hi Monocytes)",  "log(HLA-DRlo CD163hi Monocytes)", "log(HLA-DRlo S100Ahi Monocytes)", "log(Non-classical Monocytes)", "log(mDCs)", "log(pDCs)", "log(NK cells)", "log(Plasmablasts)"), each=3))



dat_plot_3g_caco <- dat_plot_3g %>% 
  filter(term != "days_post_symptoms") %>% 
  mutate(group = rep(c("mild vs control", "severe vs control"), times=9))
dat_plot_3g_caco$pos <- c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2, 6.8, 7.2, 7.8, 8.2, 8.8, 9.2)

dat_plot <- dat_plot_3g_caco %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_3g_caco$text <- res1

ggplot(data=dat_plot_3g_caco, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high), col=rev(group))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-3, 15))+
  geom_hline(yintercept=0, lty=2) +
  scale_color_manual(values=cols[2:3])+
  xlab("") +
  ylab("difference between cases and control (95% CI)") +
  theme_bw(base_size = 20)+ 
  theme(legend.position=c(0.5,0.94), legend.title=element_blank(), legend.direction = "horizontal")+
  scale_x_discrete(limits = rev(c("", unique(dat_plot_3g_caco$outcome))), breaks=unique(dat_plot_3g_caco$outcome))+
  geom_text(aes(x = pos, y = rep(13, each=nrow(dat_plot_3g_caco))), size=4.5,label = rev(res1))


dat_plot_3g_days <- dat_plot_3g %>% 
  filter(term == "days_post_symptoms")

dat_plot_3g_days$pos <- 1:nrow(dat_plot_3g_days)
dat_plot <- dat_plot_3g_days %>% dplyr::select(estimate, conf.low, conf.high) %>% 
  as.matrix(.) %>%
  formatC(., format="f", digits=2)
res1 <- paste0(dat_plot[,1], " (", dat_plot[,2], ";", dat_plot[,3], ")")
dat_plot_3g_days$text <- res1

ggplot(data=dat_plot_3g_days, aes(x=pos, y=rev(estimate), ymin=rev(conf.low), ymax=rev(conf.high))) +
  geom_pointrange(size=1.2)+
  coord_flip(clip = 'off') +
  scale_y_continuous(limits=c(-.5, 0.75))+
  geom_hline(yintercept=0, lty=2) +
  xlab("") +
  ylab("effect of days post symptom-onset (95% CI)") +
  theme_bw(base_size = 20)+ 
  scale_x_discrete(limits = rev(c(dat_plot_3g_days$outcome)), breaks=dat_plot_3g_days$outcome)+
  geom_text(aes(x = pos, y = rep(0.5, each=nrow(dat_plot_3g_days))), size=4.5,label = rev(res1))
