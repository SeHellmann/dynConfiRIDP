### Check, ob 2DSD und dynaViTE überhaupt unterschiedliche Vorhersagen machen
library(dynConfiR)
library(tidyverse)

beta_cue <- 2
beta_cue_z <- 0.3
beta_prev <- 1
beta_stim <- 2

param <- expand.grid(model=c("dynaViTE", "2DSD"),
                     effect = c("z_bias", "v_bias"),
                     stim = c(0, 1),
                     cue=c("low", "high"),
                     preveffect = c(-2,-1,1,2))
param <- param %>%
  mutate(a = 2, t0 = 0.1, sz = 0.1, sv = .2, tau = 0.7, w = ifelse(model=="2DSD", 1, 0.2), svis = 1, sigvis = 0.2, lambda = 1) %>%
  mutate(v = -4 + beta_stim*stim +  beta_cue*(cue=="high")*(effect=="v_bias") + beta_prev*preveffect,
         z = 0.5 + beta_cue_z*(ifelse(cue=="low", -1, 1))*(effect=="z_bias"))

sim_data <- param %>% group_by(model, effect, stim, cue, preveffect) %>%
  reframe(rdynaViTE(10000, a, v, t0, z, 0, sz, sv, 0, tau, w, abs(v), sigvis, svis, lambda))
table(sim_data$response)

ggplot(sim_data, aes(x=interaction(cue, preveffect), color=interaction(stim, response), y=conf))+
  geom_boxplot() +
  facet_wrap(model~effect)
ggplot(sim_data, aes(x=interaction(cue, response), color=interaction(stim), y=conf, linetype=model))+
  geom_boxplot() +
  facet_wrap(preveffect~effect)

sim_data <- sim_data %>% mutate(coherent= ifelse(ifelse(cue=="low", 0, 1)==stim, "coh", "incoh"))
distinct(sim_data[,c("cue", "stim", "coherent", "preveffect")])


agg_sim <- sim_data %>% group_by(preveffect, coherent, model, stim, response) %>%
  reframe(conf=mean(conf))
ggplot(agg_sim, aes(x=as.factor(preveffect), y=conf, color=coherent, shape=model, linetype=model,
                     group=interaction(coherent,model)))+
  geom_line()+geom_point()+
  facet_wrap(response~stim)




## Compare to data
behav = read.csv("Spielplatz/Forster_EEGWEV/prepro_behav_data_expecon1_2.csv")
# rename the variables for the model syntax (could also be changed in the model fitting function)
names(behav)[names(behav) == "sayyes"] <- "response" # 0 is no response
names(behav)[names(behav) == "respt1"] <- "rt" # cleaned, rt < 0.1 and > 2.5 s already kicked out
names(behav)[names(behav) == "conf"] <- "rating" # 0 is low confidence in being correct
names(behav)[names(behav) == "ID"] <- "sbj"
names(behav)[names(behav) == "isyes"] <- "stimulus" # 0 is catch trial (noise only)
behav$cue <- as.factor(behav$cue)

# kick out columns that contain only nan or no info at all
behav <- subset(behav, select = -c(X))
# kick out nan trials
behav = na.omit(behav)
# make sure we kicked them out
sum(is.na(behav))
head(behav)
behav %>% group_by(sbj) %>%
  reframe(nints = length(unique(intensity)))
MConf <- behav %>% group_by(sbj, study, cue, correct) %>%
  reframe(MeanConf = mean(rating))
MConf
ggplot(MConf, aes(x=as.factor(cue), shape=correct, color=as.factor(correct), group=interaction(correct, sbj), y=MeanConf))+
  geom_point() + geom_line()+
  facet_wrap(.~study)
ggplot(MConf, aes(x=interaction(cue), color=as.factor(correct),y=MeanConf))+
  geom_boxplot()+
  facet_wrap(.~study)
