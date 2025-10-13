# Source functions and clean enviroment
rm(list = ls())
source("LassoFunctions.R")
library(testthat)

# Test standardizeXY
n <- 50
p <- 3
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- rnorm(n)

out <- standardizeXY(X, Y)
Xtilde <- out$Xtilde
Ytilde <- out$Ytilde

colMeans(Xtilde)
t(Xtilde) %*% Xtilde / n
mean(Ytilde)

# Test soft
softthresh <- function(x, lambda){
  # [ToDo] Fill in to return S(x, lambda)
  return(sign(x) * max(abs(x) - lambda, 0))
} # code from class

a <- rnorm(1)
lambda <- runif(1)
softthresh(a, lambda) == soft(a, lambda)

# Test fitLASSO functions & cvLASSO function
library(hdi)
data(riboflavin)

class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

X = as.matrix(riboflavin$x)
Y = riboflavin$y

out <- standardizeXY(X, Y)
Xtilde <- out$Xtilde
Ytilde <- out$Ytilde

fit1 <- fitLASSOstandardized(Xtilde, Ytilde, lambda = 0.01)
sum(abs(fit1$beta))
fit1$fmin

fit2 <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001)
plot(fit2$fmin_vec)
sum(abs(fit2$beta_mat[, 59]))

fit3 <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = seq(from = 1, to = 0.01, by = -0.05))
plot(fit3$fmin_vec)
plot(fit3$lambda_seq)

fit4 <- fitLASSO(X, Y)
plot(fit4$lambda_seq)

fit <- cvLASSO(X, Y)
plot(fit$lambda_seq)
fit$cvm
fit$lambda_min
