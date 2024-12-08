## -----------------------------------------------------------------------------
library(SA24204172)
library(pheatmap)
n_samples <- 50
p_size <- 24
q_size <- 7
rho_value <- 0.8
sigma2 <- AR1_covariance_matrix(p_size, rho_value)
pheatmap(sigma2,cluster_rows = FALSE, cluster_cols = FALSE )

## -----------------------------------------------------------------------------
set.seed(2024)
sigma1 <- generate_symmetric_positive_definite_matrix(q_size)
pheatmap(sigma1,cluster_rows = FALSE, cluster_cols = FALSE )

## -----------------------------------------------------------------------------
set.seed(2024)
sigma <- kronecker(sigma2, sigma1)
simulated_data <- simulate_data(n_samples, p_size, q_size, sigma)

Cov <- cov(simulated_data)
pheatmap(Cov,cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 2, cellheight = 2 )

## -----------------------------------------------------------------------------
Cov_permuted <- permute(Cov,p_size,p_size, q_size, q_size)
result <- naive_estimator(Cov_permuted, p_size, q_size)
pheatmap(result$Sigma2,cluster_rows = FALSE, cluster_cols = FALSE )
pheatmap(result$Sigma1,cluster_rows = FALSE, cluster_cols = FALSE )

## -----------------------------------------------------------------------------
result_2<-propose_estimator(Cov, p_size, q_size, 0.05,50)

## -----------------------------------------------------------------------------
pheatmap(result_2$Sigma2,cluster_rows = FALSE, cluster_cols = FALSE )
pheatmap(result_2$Sigma1,cluster_rows = FALSE, cluster_cols = FALSE )

## -----------------------------------------------------------------------------
library(SA24204172)
library(pheatmap)

n_samples <- 50
p_size <- 24
q_size <- 7
rho_value <- 0.8
sigma2 <- AR1_covariance_matrix(p_size, rho_value)
pheatmap(sigma2,cluster_rows = FALSE, cluster_cols = FALSE )

