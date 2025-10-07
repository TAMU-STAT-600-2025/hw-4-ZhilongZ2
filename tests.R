# Source functions and clean enviroment
rm(list = ls())
source("LassoFunctions.R")

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

