# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi")
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
# ?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
fit1 <- fitLASSO(X = X, Y = Y)

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
nnz <- colSums(fit1$beta_mat != 0)

plot(fit1$lambda_seq, nnz, type = "b", log = "x",
     xlab = expression(lambda), ylab = "Number of nonzero coefficients",
     main = "Model size vs tuning parameter")

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark)

timing <- microbenchmark(
  fitLASSO(X, Y),
  times = 10
)
summary(timing)

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
# Median time: 1.849194 sec

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
cvfit <- cvLASSO(X = X, Y = Y, n_lambda = 30)

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
# Plot CV(lambda) vs lambda

plot(log(cvfit$lambda_seq), cvfit$cvm, type = "b", pch = 20,
     xlab = expression(log(lambda)), ylab = "CV MSE",
     main = "Cross-validation curve") 

abline(v = log(cvfit$lambda_min), col = "red", lty = 2)
abline(v = log(cvfit$lambda_1se), col = "blue", lty = 2)

legend("topleft", legend = c("lambda_min", "lambda_1se"),
       col = c("red", "blue"), lty = 2, bty = "n")
# Note: I plot with log lambda here for readability, which is reasonable due to our approach for generating lambda_seq

