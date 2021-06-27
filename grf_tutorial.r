library(grf)

# Generate data.
n <- 2000
p <- 10
# Random numbers using normal dist
# 2000 rows, 10 columns
X <- matrix(rnorm(n * p), n, p)

# 100 rows, 10 columns
X.test <- matrix(0, 101, p)

# Column 1, assign values spaced out from -2 to 2
X.test[, 1] <- seq(-2, 2, length.out = 101)

# Binary treatments
# rbinom(# obs, # trials per obs, prob success)
# 2000 obs, 1 trial per obs, .4 or .6 based on 
W <- rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))

# Outcome
# pmax/min with 0, vector level compare 
# If you look at the hist, outcomes are normally distributed but different for diff people 
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
hist(Y)

# Train a causal forest.
tau.forest <- causal_forest(X, Y, W)

# Since 1 and 3 were in Y, makes sense they're most important
tau.forest

# Estimate treatment effects for the training data using out-of-bag prediction.
# I have notes on how this works in the grf pdf
tau.hat.oob <- predict(tau.forest)
hist(tau.hat.oob$predictions)

# Estimate treatment effects for the test sample.
# This one is actually out of it
# what should the x value be?
tau.hat <- predict(tau.forest, X.test)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 2)

# Estimate the conditional average treatment effect on the full sample (CATE).
average_treatment_effect(tau.forest, target.sample = "all")

# Estimate the conditional average treatment effect on the treated sample (CATT).
average_treatment_effect(tau.forest, target.sample = "treated")

# Estimate the conditional average treatment effect on the control sample (CATC).
average_treatment_effect(tau.forest, target.sample = "control")


# Add confidence intervals for heterogeneous treatment effects; 
# growing more trees is now recommended.
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
# ?
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)

# Similar plot to the first but now there's a confidence interval ig
# why the heck 1.96?
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 1)

test_calibration(tau.forest)

# Extract the first tree from the fitted forest.
tree <- get_tree(tau.forest, 1)
# print(tree)
plot(tree)

# ?? get_leaf_node is missing
# Get a vector of node numbers for each sample.
#get_leaf_node(tree, X.test[1:10])

# Get a list of samples per node.
#get_leaf_node(tree, X.test, node.id = FALSE)
