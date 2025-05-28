# run confidence, reaction time and response based DDM
# work in progress
# author: Carina Forster, CarinaFo@github.com


# https://github.com/SeHellmann/dynConfiR

library("dynConfiR")
library("tidyverse")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library("multcomp")

## SET COLORS FOR FIGURES ####
cols = brewer.pal(5, "PuOr")
col_correct = cols[1]
col_incorrect = cols[5]
col_prevyes = cols[2]
col_prevno = cols[4]

####################################load data#####################################################
# set base directory
setwd("Spielplatz/Forster_EEGWEV/")

behav_path = file.path("prepro_behav_data_expecon1_2.csv")

# contains data from both studies
behav = read.csv(behav_path)

# rename the variables for the model syntax (could also be changed in the model fitting function)
names(behav)[names(behav) == "sayyes"] <- "response" # 0 is no response
names(behav)[names(behav) == "respt1"] <- "rt" # cleaned, rt < 0.1 and > 2.5 s already kicked out
names(behav)[names(behav) == "conf"] <- "rating" # 0 is low confidence in being correct
names(behav)[names(behav) == "ID"] <- "sbj"
names(behav)[names(behav) == "isyes"] <- "stimulus" # 0 is catch trial (noise only)

# kick out columns that contain only nan or no info at all
behav <- subset(behav, select = -c(X))

# kick out nan trials
behav = na.omit(behav)

# make sure we kicked them out
sum(is.na(behav))

# can"t fit parameters so split the data into high vs. low prob trials and fit separate models
# for each condition
high_prob = behav[behav$condition == 0.75,]
low_prob = behav[behav$condition == 0.25,]

# same for previous choice
prev_yes = behav[behav$prevresp == 1,]
prev_no = behav[behav$prevresp == 0,]

sample_data = behav[behav$sbj == sample(unique(behav$sbj), 1),]

# test DDM+Conf Model for a sample participant
start.time <- Sys.time()

# now fit the model
sample = fitRTConf(sample_data, "dynaViTE", condition="cardiac_phase",
                   opts=list(nAttempts=2, nRestarts=2),
                   fixed=list(sigvis=0, svis=1),
                   useparallel = TRUE, n.cores=4)

end.time <- Sys.time()

time.taken <- round(end.time - start.time,2)
time.taken


preds_Conf <- predictConf(sample, "dynWEV", simult_conf=TRUE)

preds_RT <- predictRT(sample, "dynWEV", subdivisions=200,
                       minrt=sample$tau+sample$t0, simult_conf = TRUE,
                      scaled=TRUE, DistConf = preds_Conf)

# plot predicitions and real data
preds_Conf$rating <- factor(preds_Conf$rating, labels=c("unsure", "sure"))
preds_RT$rating <- factor(preds_RT$rating, labels=c("unsure", "sure"))

ggplot(preds_Conf, aes(x=interaction(rating, response), y=p))+
   geom_bar(stat="identity")+
    facet_grid(cols=vars(stimulus), rows=vars(condition), labeller = "label_both")

ggplot(preds_RT, aes(x=rt, color=interaction(rating, response), y=densscaled))+
     geom_line(stat="identity")+
     facet_grid(cols=vars(stimulus), rows=vars(condition), labeller = "label_both")+
     theme(legend.position = "bottom")+ ggtitle("Scaled Densities")

ggplot(aggregate(dens~rt+correct+rating+condition, preds_RT, mean),
         aes(x=rt, color=rating, y=dens))+
     geom_line(stat="identity")+
     facet_grid(cols=vars(condition), rows=vars(correct), labeller = "label_both")+
     theme(legend.position = "bottom")+ ggtitle("Non-Scaled Densities")


# save parameters for each participant and then compare between conditions with paired t-tests
# or fit model for both conditions in one go with conditional parameters


# stimulus and response must be coded as 0 and 1
calc_sdt_params <- function(response, stimulus) {
  # Calculate hit rate, false alarm rate, and overall mean
  hit_rate <- sum(response == 1 & stimulus == 1) / sum(stimulus == 1)
  false_alarm_rate <- sum(response == 1 & stimulus == 0) / sum(stimulus == 0)
  overall_mean <- (sum(response == 1) + sum(stimulus == 1)) / length(response)

  # Apply Hautus correction to prevent undefined values
  hit_rate = ifelse(hit_rate == 1, 1 - 1 / (2 * sum(stimulus == 1)), hit_rate)
  hit_rate = ifelse(hit_rate == 0, 1 / (2 * sum(stimulus == 1)), hit_rate)
  false_alarm_rate = ifelse(false_alarm_rate == 1, 1 - 1 / (2 * sum(stimulus == 0)), false_alarm_rate)
  false_alarm_rate = ifelse(false_alarm_rate == 0, 1 / (2 * sum(stimulus == 0)), false_alarm_rate)


  # Calculate Z-scores
  z_hit <- qnorm(hit_rate)
  z_false_alarm <- qnorm(false_alarm_rate)
  z_overall_mean <- qnorm(overall_mean)

  # Calculate d' and criterion
  d_prime <- z_hit - z_false_alarm
  criterion <- -0.5 * (z_hit + z_false_alarm)

  # Return the results as a list
  return(list(d_prime = d_prime, criterion = criterion))
}

# calc_sdt_params function is at the end of the script

# Calculate d' and criterion for each participant, previous response, and cue condition
results <- behav %>%
  group_by(sbj, prevresp, cue, study) %>%
  summarize(d_prime = calc_sdt_params(response, stimulus)$d_prime,
            criterion = calc_sdt_params(response, stimulus)$criterion)

# convert to factor
results$prevresp = as.factor(results$prevresp)
results$cue = as.factor(results$cue)
results$study = as.factor(results$study)

# Plot boxplots for d' and criterion based on cue condition, previous response, and study
plot_d_prime <- ggplot(results, aes(x = cue, y = d_prime, fill = prevresp)) +
  geom_boxplot(position = position_dodge(0.8)) +
  facet_grid(~ study) +
  labs(title = "Boxplot of d' by Cue Condition, Previous Response, and Study",
       x = "Cue Condition", y = "d'") +
  scale_fill_manual(values = c(col_prevno, col_prevyes))

plot_criterion <- ggplot(results, aes(x = cue, y = criterion, fill = prevresp)) +
  geom_boxplot(position = position_dodge(0.8)) +
  facet_grid(~ study) +
  labs(title = "Boxplot of Criterion by Cue Condition, Previous Response, and Study",
       x = "Cue Condition", y = "Criterion") +
  scale_fill_manual(values = c(col_prevno, col_prevyes))

# Display the plots
plot_d_prime
plot_criterion

# split by study # 4 way interaction is total overkill and won't converge

study1 = behav[behav$study == 1,]
study2 = behav[behav$study == 2,]

prev_resp_cue_model <- glmer(response ~ cue*prevresp*stimulus+
                               (cue*prevresp*stimulus|sbj),
                             data = study1,
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)),
                             family=binomial(link='probit'))

summary(prev_resp_cue_model)


# Plot reaction time distribution for correct and error trials and high and low confidence
reaction_time_plot <- ggplot(sample_data, aes(x = rt, fill = correct)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ rating) +
  scale_fill_manual(values = c(col_correct, col_incorrect)) +
  labs(title = "Response Time Distribution for low (0) and high (1) confidence trials", x = "Response Time", y = "Density",
       fill = "Accuracy")

# Plot mean confidence based on accuracy and stimulus using a bar plot
confidence_barplot <- ggplot(sample_data, aes(x = correct, y = rating, fill = correct)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  facet_grid(. ~ stimulus) +  # Separate bars for each stimulus level
  scale_fill_manual(values = c(col_correct, col_incorrect)) +
  labs(title = "Mean Confidence Based on Accuracy and Stimulus (0 = catch trial)", x = "Accuracy", y = "Mean Confidence",
       fill = "Accuracy")

# Display the plots
reaction_time_plot
confidence_barplot

# plot mean over all participants
# plot mean reaction time for correct and error trials based on confidence level

# Calculate mean rts for each participant, accuracy, and confidence
mean_rts <- behav %>%
  group_by(sbj, correct, rating, study) %>%
  summarize(mean_rts = mean(rt))

# Plot boxplots for mean confidence based on accuracy and stimulus
rt_boxplot <- ggplot(mean_rts, aes(x = correct, y = mean_rts, fill = factor(rating))) +
  geom_boxplot() +
  facet_wrap(~ study) +
  scale_fill_manual(values = c(col_correct, col_incorrect)) +
  labs(title = "Mean response time Based on accuracy and confidence", x = "Accuracy",
       y = "Mean RTs", fill = "Confidence")

# Calculate mean confidence for each participant, accuracy, and stimulus
mean_confidence <- behav %>%
  group_by(sbj, correct, stimulus, study) %>%
  summarize(mean_rating = mean(rating))

# Plot boxplots for mean confidence based on accuracy and stimulus
confidence_boxplot <- ggplot(mean_confidence, aes(x = correct, y = mean_rating, fill = factor(stimulus))) +
  geom_boxplot() +
  facet_wrap(~ study) +
  scale_fill_manual(values = c(col_correct, col_incorrect)) +
  labs(title = "Mean Confidence Based on Accuracy and Stimulus", x = "Accuracy",
       y = "Mean Confidence", fill = "Stimulus")

# split based on
# Display the boxplot
confidence_boxplot
rt_boxplot

# fit lmer model
# transform rt to log10
behav$log10rt = log10(behav$rt)
behav$prev_conf = as.factor(behav$prevconf)
behav$congruency = as.factor(behav$congruency)

# not sure why this model does not converge
summary(lmerTest::lmer(log10rt ~ congruency +
                         (congruency|sbj), data=behav))
