#include <RcppEigen.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param Cov the number of samples
//' @param p the number of between-sample random numbers
//' @param q 111
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
List naive_estimator(Eigen::MatrixXd Cov, int p, int q) {
   // Reshape the covariance matrix into a permuted form
   Eigen::MatrixXd X_tilde = Cov; // Assume `Cov` is already correctly reshaped

   // Singular Value Decomposition (SVD)
   Eigen::JacobiSVD<Eigen::MatrixXd> svd(X_tilde, Eigen::ComputeThinU | Eigen::ComputeThinV);
   Eigen::VectorXd s = svd.singularValues();
   Eigen::MatrixXd u = svd.matrixU();
   Eigen::MatrixXd v = svd.matrixV();

   // Compute C and D
   double sqrt_s0 = std::sqrt(s(0));

   // Reshape and scale the first column of U into a p x p matrix
   Eigen::MatrixXd C = sqrt_s0 * Eigen::Map<Eigen::MatrixXd>(u.col(0).data(), p, p);

   // Reshape and scale the first row of V into a q x q matrix
   Eigen::MatrixXd D = sqrt_s0 * Eigen::Map<Eigen::MatrixXd>(v.row(0).data(), q, q);

   // Return results as a List
   return List::create(Named("Sigma2") = C, Named("Sigma1") = D);
 }

