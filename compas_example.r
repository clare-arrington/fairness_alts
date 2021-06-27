library(dplyr)
library(ggplot2)
library(lubridate)
library(grf)
library(fairness)

seed <- 42
set.seed(seed)

raw_data <- read.csv("/home/clare/SchoolWork/compas-scores-two-years.csv")
#raw_data <- read.csv("/home/clare/SchoolWork/compas-scores.csv")
nrow(raw_data)

# Charge date must be within the month
df <- dplyr::select(raw_data, age, c_charge_degree, race, age_cat, score_text, 
                    sex, juv_fel_count, juv_misd_count, juv_other_count, 
                    priors_count, days_b_screening_arrest, decile_score, 
                    is_recid, two_year_recid, c_jail_in, c_jail_out,
                    in_custody, out_custody) %>%
  filter(days_b_screening_arrest <= 30) %>%
  filter(days_b_screening_arrest >= -30) %>% na.exclude
nrow(df)

# Maybe I should just filter after?
df <- filter(df, race == 'Caucasian' | race == 'African-American')
nrow(df)

# I don't think race or sex needs to be re-leveled
df <- mutate(df, crime_factor = factor(c_charge_degree)) %>%
  mutate(age_factor = as.factor(age_cat)) %>%
  within(age_factor <- relevel(age_factor, ref = 1)) %>%
  mutate(race_factor = factor(race)) %>%
  within(race_factor <- relevel(race_factor, ref = 2)) %>%
  mutate(gender_factor = factor(sex, labels= c("Female","Male"))) %>%
  within(gender_factor <- relevel(gender_factor, ref = 2)) %>%
  mutate(score_factor = factor(score_text != "Low", labels = c("LowScore","HighScore")))

## Setup base variables for causal forest

# Y - outcome
# jail_time_duration in seconds so we convert to days
jail_time_duration <- ymd_hms(df$c_jail_out) - ymd_hms(df$c_jail_in)
#jail_time_duration <- ymd(df$out_custody) - ymd(df$in_custody)
Y <- abs(jail_time_duration) / ddays(1)
hist(Y[Y<200])

## M - merit
# 0 - no, 1 - yes
#merit <- df$two_year_recid
merit <- df$is_recid

# W - treatment
# Risk score 1 - 10
# Above 4 is high
orig_treatment <- df$score_factor

# Use this to copy same approach as first treatment for counterfactuals
get_score_text <- function(score) {
  if (score <= 4) {
    score_text <- 'LowScore'
  } else {
    score_text <- 'HighScore'
  }
  return(score_text)
}

# Affirmative action
affirm_action <- function(row) {
  if (row[['race']] == 'Caucasian') {
    risk_score <- as.numeric(row[['decile_score']]) + 1
  } else {
    risk_score <- as.numeric(row[['decile_score']]) - 1
  }
  return(risk_score)
}

# Sapply doesn't need to specify level; good for vector
affirm_scores <- apply(df, MARGIN = 1, FUN = affirm_action)
affirm_text <- sapply(affirm_scores, FUN = get_score_text)
affirm_treatment <- relevel(factor(affirm_text), ref = 2)

# Perfect prediction
# Add 1 just to match how it was factored in the first case
# Not in factor mode but that should be fine?
perfect_treatment <- df$two_year_recid + 1

## Setup features
X <- df[,c("juv_fel_count", "juv_misd_count", "juv_other_count", 
           "priors_count")] 

# For the factors, I have to convert to numeric first.
# I guess I'll have to ref back to the df if needed
X$race <- as.numeric(df$race_factor) - 1
X$gender <- as.numeric(df$gender_factor) - 1
X$age <- as.numeric(df$age_factor) - 1
X$crime <- as.numeric(df$crime_factor) 

# Sample data.
dim(X)
num_samples <- 5000
train_sample = sample(1:5278, num_samples, replace=FALSE)

X.train <- X[train_sample,]
X.test <- X[-train_sample,]

#X.train <- as.matrix(X[train_sample,]$race)
#X.test <- as.matrix(X[-train_sample,]$race)

#length(Y)
Y.train <- Y[train_sample]
Y.test <- Y[-train_sample]

# Change this for different counter facts
W.train <- as.numeric(orig_treatment)[train_sample] - 1
W.test <- as.numeric(orig_treatment)[-train_sample] - 1

# Change this for different attributes
A.train <- X$race
A.test <- X$race

# Before assigning the data to causal_forest(), Athey et al. (2019) 
# recommend to fit separate regression trees to estimate Y and W hat
# These are estimates of m(X) = E[Y | X]
forest.Y <- regression_forest(X=X.train, Y=Y.train, seed=seed)
Y.hat <- forest.Y$predictions
hist(Y.hat)

# These are estimates of the propensity score E[W | X]
W.hat <- regression_forest(X=X.train, Y=W.train, seed=seed)$predictions
hist(W.hat)

# Compare variable importance
forest.Y.varimp <- variable_importance(forest.Y)
selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)
print('Selected Variables are:')
selected.vars
#forest.Y.varimp

#dimnames(X)[2]
# Order of importance for orig: priors count, race, age, gender, crime, juv_other

# Train a causal forest.
tau.forest <- causal_forest(X=X.train[, selected.vars], Y=Y.train, W=W.train,
                            Y.hat=Y.hat, W.hat=W.hat, seed = seed)

# Asses a forest’s goodness of fit. 
# 1 for mean suggests that the mean forest prediction is correct 
# 1 for differential suggests that heterogeneity is captured in the underlying signal
# test_calibration(tau.forest)

# See where the splitting happens
# A matrix of split depth by feature index, 
# where each value is the number of times the feature was split on at that depth.
print('Split frequencies are:')
split_frequencies(tau.forest)

# Estimate the conditional average treatment effect on the types of samples.
average_treatment_effect(tau.forest, target.sample = "all") # CATE
average_treatment_effect(tau.forest, target.sample = "treated") # CATT
print("Control ATE")
average_treatment_effect(tau.forest, target.sample = "control") # CATC
print("Overlap ATE")
average_treatment_effect(tau.forest, target.sample = "overlap") # ??

# CATE: tau(x) = E[Y(1) - Y(0) | X = x]
# Estimate treatment effects for the training data using out-of-bag prediction.
tau.hat.train <- predict(tau.forest)$predictions 
hist(tau.hat.train, 
     main = 'Outcome Predictions for Out-of-Bag Training Samples',
     xlab = 'Estimated Treatment Effect')

X.train$treatment.effect <- tau.hat.train

sprintf("Mean ETE for White Defendants = %f", mean(subset(X.train, race == 0)$treatment.effect))
sprintf("Mean ETE for Black Defendants = %f", mean(subset(X.train, race == 1)$treatment.effect))

sprintf("Mean ETE for Both = %f", mean(X.train$treatment.effect))
sprintf("IQR ETE for Both = %f", IQR(X.train$treatment.effect))
sprintf("SD ETE for Both = %f", sd(X.train$treatment.effect))

X.train$Ethnicity <- ifelse(X.train$race == 0, 'Caucasian', 'African-American')

ggplot(X.train, aes(x=treatment.effect, fill = Ethnicity, color = Ethnicity)) + 
  geom_histogram(position="identity", alpha=.3, binwidth = 5) + 
  scale_color_manual(aesthetics=c("colour", "fill"), 
                     values=c(rgb(.8,.1,.1), rgb(.1,.1,.8))) +
  xlab('Estimated Treatment Effect') + ylab('Frequency') + 
  ggtitle('Outcome Predictions for Out-of-Bag Training Samples') +
  theme(text = element_text(size = 15)) 

# E[Y | X, W = 0]
mu.hat.train.0 <- (Y.hat - W.hat * tau.hat.train)
# E[Y | X, W = 1]
mu.hat.train.1 <- (Y.hat + (1 - W.hat) * tau.hat.train)

mu.hat <- data.frame( jail.time = c(mu.hat.train.0, mu.hat.train.1), 
                      Treatment = c(rep("Won't Reoffend", num_samples), 
                                    rep('Will Reoffend', num_samples)))

ggplot(mu.hat, aes(x=jail.time, fill = Treatment, color = Treatment)) + 
  geom_histogram(position="identity", alpha=.3, binwidth = 5) + 
  scale_color_manual(aesthetics=c("colour", "fill"), 
                     values=c(rgb(.8,.1,.1), rgb(.1,.1,.8))) +
  xlab('Jail Time (Days)') + ylab('Frequency') + 
  ggtitle('Outcome Predictions for Out-of-Bag Training Samples') +
  theme(text = element_text(size = 15)) 

# Estimate treatment effects for the test sample.
# This one is actually out of it
tau.hat.test <- predict(tau.forest, X.test[, selected.vars])
hist(tau.hat.test$predictions, 
     main = 'Outcome Predictions for Testing Samples',
     xlab = 'Estimated Treatment Effect')

ggplot(mu.hat, aes(x=jail.time, fill = Treatment, color = Treatment)) + 
  geom_histogram(position="identity", alpha=.3) + 
  scale_color_manual(aesthetics=c("colour", "fill"), 
                     values=c(rgb(1,0,0), rgb(0,0,1))) +
  xlab('Jail Time (Days)') + ylab('Frequency') + 
  ggtitle('Outcome Predictions for Out-of-Bag Training Samples') +
  theme(text = element_text(size = 15)) 

# E[Y | X, W = 0]
#mu.hat.test.0 <- (Y.hat - W.hat * tau.hat.test)$predictions
#hist(mu.hat.test.0)

# E[Y | X, W = 1]
#mu.hat.test.1 <- (Y.hat + (1 - W.hat) * tau.hat.test)$predictions
#hist(mu.hat.test.1)

# use d eps nu

# to get n(x)
# start with baseline dist
# focus on nu from that

e.hat <- W.hat

n <- 8
p <- num_samples

IPW <- ifelse(W.train == 1, 1 / e.hat, 1 / (1 - e.hat))

plot.df <- data.frame(value = as.vector(X.train),
                      variable = colnames(X.train)[rep(1:p, each = n)],
                      W = as.factor(W.train),
                      IPW = IPW)

ggplot(plot.df, aes(x = value, weight = IPW, fill = W.train)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  facet_wrap( ~ variable, ncol = 2)


## Get Fairness
# X.df <- as.data.frame(X)
# X.df$outcome <- merit
# X.df$treatment <- as.numeric(affirm_treatment) - 1
# 
# preds <- pred_rate_parity(X.df, outcome='outcome', preds='treatment',
#                           group='race', base='0')
# 
# preds$Metric

## Check data
# black <- subset(X.df, race == 1)
# (subset(black, treatment == 1)) / nrow(black)
# black <- subset(df, race_factor == 'African-American')
# 
# temp <- ymd_hms(black$c_jail_out) - ymd_hms(black$c_jail_in)
# mean(temp) / ddays(1)
# white <- subset(X.df, race == 0)
# nrow(subset(white, treatment == 1)) / nrow(white)
# white <- subset(df, race_factor == 'Caucasian')
# 
# temp <- ymd_hms(white$c_jail_out) - ymd_hms(white$c_jail_in)
# mean(temp) / ddays(1)
