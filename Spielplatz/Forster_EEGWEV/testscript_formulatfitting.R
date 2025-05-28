devtools::document()
devtools::load_all()
devtools::install()
library(dynConfiRIDP)
rm(list=ls())
setwd("../..")
# contains data from both studies
behav = read.csv("Spielplatz/Forster_EEGWEV/prepro_behav_data_expecon1_2.csv")

# rename the variables for the model syntax (could also be changed in the model fitting function)
names(behav)[names(behav) == "sayyes"] <- "response" # 0 is no response
names(behav)[names(behav) == "respt1"] <- "rt" # cleaned, rt < 0.1 and > 2.5 s already kicked out
names(behav)[names(behav) == "conf"] <- "rating" # 0 is low confidence in being correct
names(behav)[names(behav) == "ID"] <- "sbj"
names(behav)[names(behav) == "isyes"] <- "stimulus" # 0 is catch trial (noise only)
behav$cue <- as.factor(behav$cue)
behav$X
# kick out columns that contain only nan or no info at all
behav <- subset(behav, select = -c(X))
# kick out nan trials
behav = na.omit(behav)
# make sure we kicked them out
sum(is.na(behav))
table(behav$sbj, behav$study)

head(behav)
# manipulations = list(v~intensity+cue, z~cue)
# manipulations = list(v~intensity, d~cue, z~cue)
# manipulations = list(v~stimulus+cue+prevresp*prevconf)
# fixed = list(sym_thetas = FALSE, s=1, svis=1,sigvis=0, d=0)

# model = "dynaViTE";nRatings = 2; restr_tau =Inf; precision=2;logging=TRUE; opts=list(); optim_method = "bobyqa";
# useparallel = TRUE; n.cores=5

setwd("Spielplatz/Forster_EEGWEV/")
data <- behav[behav$study==1, c("sbj", "stimulus", "cue", "intensity", "response", "rating", "rt", "prevresp", "prevconf")]
data <- subset(data, sbj==22)
res_22 <- fitRTConf_formula(data, "dynaViTE",
                            manipulations=list(v~stimulus+cue+prevresp*prevconf),
                            fixed= list(sym_thetas = FALSE, s=1, svis=1,sigvis=0, d=0),
                            nRatings=2, restr_tau = Inf, precision=2, logging=TRUE,
                            opts=list(maxit=2000),
                            optim_method = "Nelder-Mead", useparallel = TRUE,
                            n.cores = 5)
save(res_22, file="first_fit_res_sbj22.RData")



# manipulations = list(v~intensity+cue, z~cue)
# manipulations = list(v~intensity, d~cue, z~cue)
# manipulations = list(v~stimulus+cue+prevresp*prevconf)
# fixed = list(sym_thetas = FALSE, s=1, svis=1,sigvis=0, d=0)

# model = "dynaViTE";nRatings = 2; restr_tau =Inf; precision=2;logging=TRUE; opts=list(); optim_method = "bobyqa";
# useparallel = TRUE; n.cores=5


setwd("Spielplatz/Forster_EEGWEV/")

load("preliminary_fit.RData")

# Or:
{

  data <- behav[behav$study==1, c("sbj", "stimulus", "cue", "intensity", "response", "rating", "rt", "prevresp", "prevconf")]
  unique(data$sbj)
  data <- subset(data, sbj%in% c(10,11))
  res_1011 <- fitRTConfModels_formula(data, "dynaViTE", manipulations, fixed=fixed,
                                      nRatings=2, restr_tau = Inf, precision=2, logging=TRUE,
                                      opts=list(maxit=2000, nAttempts=2, nRestarts=3),
                                      optim_method = "Nelder-Mead",
                                      parallel="both", n.cores=c(2,2))
  save(res_1011, file="first_fit_res_sbj1011.RData")




  setwd("Spielplatz/Forster_EEGWEV/")
  load("first_fit_res_sbj1011.RData")
  load("first_fit_res_sbj22.RData")
  load("first_fit_res_sbj7.RData")
  # load("first_fit_res.RData") # same as ..._sbj7
  res_1011[[1]]$sbj
  res_1011[[1]]$beta
  res_1011[[2]]$beta
  res_22$beta
  res_7 <- res
  res_7$beta
  res_7$model_matrix

  simu <- simulate_dynaViTE_formula(res_7, n=10000)
  simu

  # data <- behav[behav$study==1, c("sbj", "stimulus", "cue", "intensity", "response", "rating", "rt", "prevresp", "prevconf")]
  # data <- subset(data, sbj==7)
  # res_7$maxt0 <- min(data$rt)
  # res_7$restr_tau <- Inf
  #
  # beta <- res_7
  #
  # pred_matrix=NULL; fixed=NULL; #method="simulation_matrix";
  # model_matrix=NULL; n=NULL; maxt0=NULL; restr_tau=NULL;
  # delta=0.01; maxrt=15; simult_conf = FALSE;
  # process_results=FALSE; seed=NULL
  # n <- 10000


  res_22$sbj <- 22
  res_22$model <- "dynaViTE"
  res_7$sbj <- 7
  res_7$model <- "dynaViTE"
  res_7$restr_tau <- Inf
  res_22$restr_tau <- Inf

  res_22$maxt0 <- behav %>% filter(sbj==22) %>% select(rt) %>% min(.)
  res_22$maxt0 <- behav %>% filter(sbj==22) %>% select(rt) %>% min(.)

  res_1011[[1]]$restr_tau <- Inf
  res_1011[[2]]$restr_tau <- Inf
  res_1011[[1]]$maxt0 <- behav %>% filter(sbj==res_1011[[1]]$sbj) %>% select(rt) %>% min(.)
  res_1011[[2]]$maxt0 <- behav %>% filter(sbj==res_1011[[2]]$sbj) %>% select(rt) %>% min(.)

  res <- c(res_1011, list(res_22), list(res_7))

  simus <- lapply(res, simulate_dynaViTE_formula, n=10000)


  structure(simus[[1]])
  do.call(rbind, simus)
  names(res_1011[[1]])
  names(res_22)

  predictions <- data.frame()
  i <- 1
  for (i in 1:length(simus)) {
    temp <- as.data.frame(simus[[i]], row.names = 1:nrow(simus[[i]]))
    #names(temp)
    temp$sbj <- res[[i]]$sbj
    if (nrow(predictions)>0) {
      predictions[,setdiff(names(temp), names(predictions))] <- NA
      temp[,setdiff(names(predictions), names(temp))] <- NA
    }
    predictions <- rbind(predictions, temp)
  }

}

  head(predictions)
library(tidyverse)
agg_predictions <- predictions %>%
  mutate(response = response*rating) %>%
  group_by(sbj, stimulus, cue0.75, prevconf, prevresp) %>%
  mutate(N=n()) %>%
  group_by(sbj, stimulus, cue0.75, prevconf, prevresp, response) %>%
  reframe(p = n()/N[1])
agg_predictions <- filter(agg_predictions, is.na(stimulus)) %>% select(-stimulus) %>%
  merge(data.frame(stimulus=c(0,1))) %>%
  rbind(.,filter(agg_predictions, !is.na(stimulus)))

agg_data <- behav %>% filter(sbj %in% c(7, 22, 10, 11)) %>%
  mutate(response= ifelse(response==1, 1, -1)*(rating+1),
         cue0.75 = ifelse(cue==0.75, 1, 0)) %>%
  group_by(sbj, stimulus, cue0.75, prevconf, prevresp) %>%
  mutate(N=n()) %>%
  group_by(sbj, stimulus, cue0.75, prevconf, prevresp, response) %>%
  reframe(p = n()/N[1])
ggplot(agg_predictions, aes(x=interaction(stimulus, as.factor(response)), y=p))+
  geom_bar(data = agg_data, stat="identity")+
  geom_point()+
  facet_grid(sbj~cue0.75+prevresp+prevconf)
agg_data <- agg_data %>% mutate(correct=as.numeric(ifelse(response<0, 0, 1)==stimulus))
agg_predictions <- agg_predictions %>% mutate(correct=as.numeric(ifelse(response<0, 0, 1)==stimulus))
pdf("data_predictions_per_sbj.pdf")
for (i in unique(agg_predictions$sbj)) {
  p <- ggplot(subset(agg_predictions,sbj==i), aes(x=as.factor(response), y=p, fill=as.factor(correct)))+
    geom_bar(data = subset(agg_data, sbj==i),stat="identity")+
    geom_point()+ggtitle(i)+
    facet_grid(prevresp+prevconf~cue0.75+stimulus, labeller=label_both)
  show(p)
}
dev.off()







agg_predictions <- predictions %>%
  mutate(response = response*rating) %>%
  group_by(sbj, stimulus, cue0.75) %>%
  mutate(N=n()) %>%
  group_by(sbj, stimulus, cue0.75, response) %>%
  reframe(p = n()/N[1])
agg_predictions <- filter(agg_predictions, is.na(stimulus)) %>% select(-stimulus) %>%
  merge(data.frame(stimulus=c(0,1))) %>%
  rbind(.,filter(agg_predictions, !is.na(stimulus)))
agg_data <- behav %>% filter(sbj %in% c(7, 22, 10, 11)) %>%
  mutate(response= ifelse(response==1, 1, -1)*(rating+1),
         cue0.75 = ifelse(cue==0.75, 1, 0)) %>%
  group_by(sbj, stimulus, cue0.75) %>%
  mutate(N=n()) %>%
  group_by(sbj, stimulus, cue0.75, response) %>%
  reframe(p = n()/N[1])
ggplot(agg_predictions, aes(x=as.factor(response), y=p, shape="Prediction"))+
  geom_bar(aes(fill="Data"),data = agg_data, stat="identity")+
  geom_point()+ xlab("Response (-1: no signal, 1: signal)  x Rating (1: unsure, 2:sure)")+
  facet_grid(sbj~stimulus+cue0.75, labeller = label_both)









save(res, predictions, simus, file="preliminary_fit.RData")
res[[4]]$beta



