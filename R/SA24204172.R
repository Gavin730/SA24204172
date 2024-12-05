#' @import MASS
#' @import boot
#' @import DAAG
#' @import bootstrap
#' @import ggplot2
#' @import pheatmap
#' @import stats
#' @import knitr
#' @import microbenchmark
#' @import rmarkdown
#' @importFrom Matrix kronecker
#' @importFrom stats runif
#' @useDynLib SA24204172
#' @importFrom Rcpp sourceCpp
NULL

#' @title A correlation matrix generator using R
#' @description generate a AR correlation matrix
#' @param size the number of variance
#' @param rho the correlation parameter
#' @return A correlation matrix
#' @examples
#' \dontrun{
#'     Cov <- AR1_covariance_matrix(size,rho)
#'     pheatmap(Cov,cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 2, cellheight = 2 )
#' }
#' @export
AR1_covariance_matrix <- function(size, rho) {
  sigma <- matrix(0, nrow = size, ncol = size)
  for (l in 1:size) {
    for (m in 1:size) {
      sigma[l, m] <- rho^abs(l - m)
    }
  }
  return(sigma)
}

#' @title A simulated data generator using R
#' @description generate  simulated data follow multivariate normal distribution
#' @param n the number of simulated data
#' @param p the row number
#' @param q the column number
#' @param sigma covariance matrix
#' @return A simulated data matrix
#' @examples
#' \dontrun{
#'     simulated_data <- simulate_data(n_samples, p_size, q_size, sigma)
#' }
#' @export
simulate_data <- function(n, p, q, sigma) {
  samples <- mvrnorm(n = n, mu = rep(0, p * q), Sigma = sigma)
  return(samples)
}

#' @title A banding operator using R
#' @description band a matrix
#' @param X matrix
#' @param k bandwidth
#' @return A banded matrix
#' @examples
#' \dontrun{
#'     data_banded <- banded(data,k)
#' }
#' @export
banded <- function(X, k) {
  for (i in 1:nrow(X)) {
    for (j in 1:ncol(X)) {
      if (abs(j - i) > k) {
        X[i, j] <- 0
      }
    }
  }
  return(X)
}

#' @title A reshape operator using R
#' @description reshapes a matrix's blocks
#' @param mat matrix
#' @param A_dim1 dim1
#' @param A_dim2 dim2
#' @param B_dim1 dim3
#' @param B_dim2 dim4
#' @return A reshaped matrix
#' @examples
#' \dontrun{
#'     Cov_permuted <- permute(Cov,k)
#' }
#' @export
permute <- function(mat, A_dim1, A_dim2, B_dim1, B_dim2) {
  ans <- matrix(0, nrow = A_dim1 * A_dim2, ncol = B_dim1 * B_dim2)
  for (j in 1:A_dim2) {
    for (i in 1:A_dim1) {
      submat <- mat[((i - 1) * B_dim1 + 1):(i * B_dim1), ((j - 1) * B_dim2 + 1):(j * B_dim2)]
      ans[(A_dim1 * (j - 1) + i), ] <- as.vector(t(submat))
    }
  }
  return(ans)
}
#' @title A symmetric positive definite matrix generator using R
#' @description generate a random symmetric positive definite matrix
#' @param p dim
#' @return A symmetric positive definite matrix
#' @examples
#' \dontrun{
#'     Cov_permuted <- permute(Cov,k)
#' }
#' @export
generate_symmetric_positive_definite_matrix <- function(p){
  A <- matrix(runif(p * p), nrow = p, ncol = p)
  return(A%*%t(A))
}

objective_function <- function(params, p, q, empirical_cov,penalty_weight) {
  # Reshape params to form Sigma1 and Sigma2
  sigma1 <- matrix(params[1:(p*p)], nrow = p, ncol = p)
  sigma2 <- matrix(params[(p*p + 1):(p*p + q*q)], nrow = q, ncol = q)

  # Kronecker product of Sigma1 and Sigma2
  cov_estimate <- kronecker(sigma1, sigma2)

  # Creating matrix D and B
  D <- diag(rep(1, p - 1))
  D <- rbind(D, rep(0, p - 1))
  B <- matrix(0, nrow = p, ncol = p - 1)

  B[1, ] <- D[1, ]
  for (i in 2:nrow(D)) {
    B[i, ] <- D[i, ] - D[i - 1, ]
  }
  D <- B

  # Calculate penalty term
  penalty_term <- sum(abs(diff(diff(sigma1, axis = 1), axis = 1)))

  # Frobenius norm and objective value
  frobenius_norm <- sqrt(sum((empirical_cov - cov_estimate)^2))
  objective_value <- frobenius_norm + penalty_weight * penalty_term

  return(objective_value)
}

#' @title proposed estimator using R
#' @description solve covariance matrix
#' @param Cov covariance matrix
#' @param p dim1
#' @param q dim2
#' @param lambda penalty weight
#' @param maxiter run times
#' @return solution
#' @examples
#' \dontrun{
#'    result <- propose_estimator(Cov,p,q,maxiter)
#' }
#' @export
propose_estimator <- function(Cov, p, q,lambda, maxiter) {
  # Initial guess (diagonal identity matrices)
  initial_guess <- c(diag(1, p), diag(1, q))

  # Use optim function to minimize with correct passing of additional parameters
  result <- optim(initial_guess,
                  fn = function(params) objective_function(params, p, q, Cov,lambda),  # Anonymous function to pass extra args
                  method = "BFGS",
                  control = list(maxit = maxiter))

  optimal_params <- result$par
  optimal_sigma1 <- matrix(optimal_params[1:(p*p)], nrow = p, ncol = p)
  optimal_sigma2 <- matrix(optimal_params[(p*p + 1):(p*p + q*q)], nrow = q, ncol = q)

  return(list(Sigma2=optimal_sigma1, Sigma1=optimal_sigma2))
}

