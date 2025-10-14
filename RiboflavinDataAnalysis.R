# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
full_fit <- fitLASSO(X, Y, n_lambda = 60, eps = 0.001)

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
num_nonzero <- colSums(full_fit$beta_mat != 0)
plot(log(full_fit$lambda_seq), num_nonzero,
     type = "l", col = "blue", lwd = 2,
     xlab = "Log(Lambda)",
     ylab = "Number of Non-Zero Coefficients",
     main = "LASSO Coefficient Path for Riboflavin Data")
grid()

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library("microbenchmark")
timing_results <- microbenchmark(
  fitLASSO(X, Y, n_lambda = 60, eps = 0.001),
  times = 10
)
print(timing_results)

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
# Median: 4.8392 seconds
# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
cv_fit <- cvLASSO(X, Y, n_lambda = 30, k = 5, eps = 0.001)

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
log_lambdas <- log(cv_fit$lambda_seq)
cvm <- cv_fit$cvm
cvse <- cv_fit$cvse

# Create the plot
plot(log_lambdas, cvm,
     ylim = range(cvm - cvse, cvm + cvse), # Ensure error bars fit
     pch = 19, col = "red",
     xlab = "Log(Lambda)",
     ylab = "Cross-Validation MSE",
     main = "5-Fold CV Error for Riboflavin Data")
grid()