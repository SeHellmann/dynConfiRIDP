##### Script to check parameter recovery for 2DSD and dynWEV

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/SeqSamplingConfidence")
# # or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

install.packages("dynConfiR_0.0.4.tar.gz", repos = NULL, type = "source")
{
  library(dynConfiR)
  library(tidyverse)
  library(parallel)
}
#=====================


behav = read.csv("prepro_behav_data_expecon1_2_anon.csv")
names(behav)
max(behav$ID)
# rename the variables for the model syntax (could also be changed in the model fitting function)
data <- behav %>% select(
  sbj = ID, study,
  stimulus = isyes,   # 0 is catch trial (noise only)
  cue,
  response = sayyes, # 0 is no response
  rt = respt1,       # cleaned, rt < 0.1 and > 2.5 s already kicked out
  rating = conf,     # 0 is low confidence in being correct
  prevresp, prevconf) %>%
  mutate(cue = as.factor(cue),
         sbj = sbj+ study*10^5) %>%
  na.omit(.) # kick out nan trials

# make sure we kicked them out
sum(is.na(data))

table(data$sbj, data$study)
head(data)
table(data$cue, data$study)
table(data$rating, data$study)
table(data$stimulus, data$study)

unique(data$sbj)

data %>% group_by(sbj) %>%
  summarise(nratings = length(unique(rating))) %>%
  filter(nratings < 2)

# manipulations = list(v~intensity+cue, z~cue)
# manipulations = list(v~intensity, d~cue, z~cue)
manipulations_driftrate = list(v~stimulus+cue)
fixed = list(sym_thetas = FALSE, s=1, svis=1, d=0)
#parallel::detectCores()
fits_driftrate <- fitRTConfModels_formula(data, "dynaViTE", manipulations_driftrate, fixed=fixed,
                                          nRatings=2, restr_tau = Inf, precision=2, logging=TRUE,
                                          opts=list(maxit=4000, nAttempts=3, nRestarts=4),
                                          optim_method = "Nelder-Mead", parallel = "both",
                                          n.cores = c(7, 3)) # 3 cores per sbj; 7 sbjs in parallel
#parallel::detectCores()
fits_driftrate <- fitRTConfModels_formula(data, "dynaViTE", manipulations_driftrate, fixed=fixed,
                                          nRatings=2, restr_tau = Inf, precision=2, logging=TRUE,
                                          opts=list(maxit=4000, nAttempts=3, nRestarts=4),
                                          optim_method = "Nelder-Mead", parallel = FALSE,
                                          n.cores = c(7, 3)) # 3 cores per sbj; 7 sbjs in parallel
save(fits_driftrate, file="fits_driftrate.RData")


