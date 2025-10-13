# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  # [ToDo] Center and scale X
  n <- nrow(X)
  p <- ncol(X)
  
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, n, p, byrow = TRUE)
  
  normsX <- colSums(Xcentered ^ 2) / n
  weights <- sqrt(normsX)
  Xtilde <- Xcentered %*% diag(1 / weights)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  sign(a) * pmax(abs(a) - lambda, 0)
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  r <- as.numeric(Ytilde - Xtilde %*% beta)
  loss <- sum(r^2) / (2 * n)
  pen  <- lambda * sum(abs(beta))
  return(loss + pen)
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde))
    stop("nrow(Xtilde) must equal length(Ytilde).")
  
  #[ToDo]  Check that lambda is non-negative
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0)
    stop("lambda must be a single non-negative number.")
  
  # Parameter set up
  Ytilde <- as.numeric(Ytilde)
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)) {
    beta <- numeric(p)
  } else {
    if (length(beta_start) != p) {
      stop("beta_start length must equal p.")
    }
    beta <- as.numeric(beta_start)
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  # Initial objective
  r <- as.numeric(Ytilde - Xtilde %*% beta)
  f_prev <- lasso(Xtilde, Ytilde, beta, lambda)
  
  max_iter <- 10000L
  for (k in 1:max_iter) {
    r <- as.numeric(Ytilde - Xtilde %*% beta)
    
    for (j in 1:p) {
      xj <- Xtilde[, j]
      z <- beta[j] + sum(xj * r) / n
      beta_new_j <- soft(z, lambda)
      r <- r + xj * (beta[j] - beta_new_j)
      beta[j] <- beta_new_j
    }
    
    f_curr <- sum(r^2) / (2 * n) + lambda * sum(abs(beta))
    if ((f_prev - f_curr) < eps){
      # Return 
      # beta - the solution (a vector)
      # fmin - optimal function value (value of objective at beta, scalar)
      return(list(beta = beta, fmin = f_curr))
    }
    f_prev <- f_curr
  }
  
  warning("Reached max_iter without meeting eps criterion; returning last iterate.")
  f_min <- sum(r^2) / (2 * n) + lambda * sum(abs(beta))
  return(list(beta = beta, fmin = f_min))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde))
    stop("nrow(Xtilde) must equal length(Ytilde).")
   
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    if (!is.numeric(lambda_seq)) {
      warning("Non-numeric lambda_seq supplied; will generate a default sequence.")
      lambda_seq <- NULL
    } else {
      lambda_seq <- unique(lambda_seq[lambda_seq >= 0])
      if (length(lambda_seq) == 0L) {
        warning("Negative lambdas supplied; will generate a default sequence.")
        lambda_seq <- NULL
      } else {
        lambda_seq <- sort(lambda_seq, decreasing = TRUE)
      }
    }
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  lambdaMax <- function(Xtilde, Ytilde){
    n <- nrow(Xtilde)
    cors <- crossprod(Xtilde, Ytilde) / n
    return(max(abs(cors)))
  }
  
  if (is.null(lambda_seq)) {
    lambda_max <- lambdaMax(Xtilde, Ytilde)
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  } else {
    n_lambda <- length(lambda_seq)
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  p <- ncol(Xtilde)
  beta_mat <- matrix(nrow = p, ncol = n_lambda)
  fmin_vec <- numeric(n_lambda)
  
  beta_start <- numeric(p)
  for (i in seq_len(n_lambda)) {
    lam <- lambda_seq[i]
    fit <- fitLASSOstandardized(Xtilde, Ytilde, lambda = lam, beta_start = beta_start, eps = eps)
    beta_mat[, i] <- fit$beta
    fmin_vec[i] <- fit$fmin
    beta_start <- fit$beta
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  std <- standardizeXY(X, Y)
  Xtilde  <- std$Xtilde
  Ytilde  <- std$Ytilde
  Ymean   <- std$Ymean
  Xmeans  <- std$Xmeans
  weights <- std$weights

  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fit <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, 
                                  n_lambda = n_lambda, eps = eps)
 
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  beta_mat <- sweep(fit$beta_mat, 1, weights, FUN = "/")
  beta0_vec <- as.numeric(Ymean - crossprod(Xmeans, beta_mat))
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = fit$lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  n <- nrow(X)
  p <- ncol(X)
  
  fit <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq_used <- fit$lambda_seq
  L <- length(lambda_seq_used)
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)) {
    fold_ids <- sample(rep_len(1:k, n))
  } else {
    if (length(fold_ids) != n) stop("fold_ids must have length n.")
    if (!identical(sort(unique(fold_ids)), 1:max(fold_ids))) stop("fold_ids should be positive integers 1..K.")
    k <- max(fold_ids)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  mse_folds <- matrix(nrow = k, ncol = L)
  
  for (fold in 1:k) {
    X_tr <- X[fold_ids != fold, , drop = FALSE]
    Y_tr <- Y[fold_ids != fold]

    X_va <- X[fold_ids == fold, , drop = FALSE]
    Y_va <- Y[fold_ids == fold]
    
    fit_tr <- fitLASSO(X_tr, Y_tr, lambda_seq = lambda_seq_used, n_lambda = L, eps = eps)
    yhat_mat <- matrix(fit_tr$beta0_vec, nrow = nrow(X_va), ncol = L, byrow = TRUE) +
      X_va %*% fit_tr$beta_mat

    resid_mat <- yhat_mat - matrix(Y_va, nrow = length(Y_va), ncol = L, byrow = FALSE)
    mse_folds[fold, ] <- colMeans(resid_mat^2)
  }
  
  cvm  <- colMeans(mse_folds)
  cvse <- apply(mse_folds, 2, sd) / sqrt(k)
  
  # [ToDo] Find lambda_min
  idx_min <- which.min(cvm)
  lambda_min <- lambda_seq_used[idx_min]

  # [ToDo] Find lambda_1SE
  thresh <- cvm[idx_min] + cvse[idx_min]
  # Largest λ with CV(λ) <= threshold
  idx_1se <- tail(which(cvm <= thresh), 1)
  lambda_1se <- lambda_seq_used[idx_1se]
  
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq_used, beta_mat = fit$beta_mat, beta0_vec = fit$beta0_vec, 
              fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, 
              cvm = cvm, cvse = cvse))
}

