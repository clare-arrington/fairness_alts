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
jail_time_duration <- ymd(df$out_custody) - ymd(df$in_custody)

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

# Sample data.
dim(X)
num_samples <- 5000
train_sample = sample(1:5278, num_samples, replace=FALSE)

X.train <- as.matrix(X[train_sample,]$race)
X.test <- as.matrix(X[-train_sample,]$race)

#length(Y)
Y.train <- Y[train_sample]
Y.test <- Y[-train_sample]

# Change this for different counter facts
W.train <- as.numeric(orig_treatment)[train_sample] - 1
W.test <- as.numeric(orig_treatment)[-train_sample] - 1

# Change this for different attributes
A.train <- X$race
A.test <- X$race

forest.Y <- regression_forest(X=X.train, Y=Y.train, seed=seed)
Y.hat <- forest.Y$predictions
hist(Y.hat)

# These are estimates of the propensity score E[W | X]
W.hat <- regression_forest(X=X.train, Y=W.train, seed=seed)$predictions
hist(W.hat)

# Train a causal forest.
tau.forest <- causal_forest(X=X.train, Y=Y.train, W=W.train,
                            Y.hat=Y.hat, W.hat=W.hat, seed = seed)

test_calibration(tau.forest)

average_treatment_effect(tau.forest, target.sample = "all") # CATE
average_treatment_effect(tau.forest, target.sample = "treated") # CATT
average_treatment_effect(tau.forest, target.sample = "control") # CATC
average_treatment_effect(tau.forest, target.sample = "overlap") # ??

tau.hat.train <- predict(tau.forest)$predictions 
hist(tau.hat.train, 
     main = 'Outcome Predictions for Out-of-Bag Training Samples',
     xlab = 'Estimated Treatment Effect')

