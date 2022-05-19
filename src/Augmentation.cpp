// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------

//' Calculate Augmentation Components
//'
//' @param covars One row per subject.
//' @param mu Grand mean.
//' @param psi Influence functions.
//' @return List.
//' @noRd 
// [[Rcpp::export]]

SEXP CalcAugComp(
	const arma::mat covars,
	const arma::rowvec mu,
	const arma::colvec psi
){
	// Subjects.
	const int n = covars.n_rows;

	// Residual matrix.
	const arma::mat resid = covars.each_row() - mu;

	// Mean resid.
	const arma::colvec mean_resid = arma::mean(resid, 0).t();

	// Gamma.
	const arma::colvec gamma = resid.t() * psi / (n * n);

	// Sigma.
	const arma::mat sigma = resid.t() * resid / (n * n);

  // Output.
	return Rcpp::List::create(
		Rcpp::Named("gamma")=gamma,
		Rcpp::Named("mean_resid")=mean_resid,
		Rcpp::Named("n")=n,
		Rcpp::Named("resid")=resid,
		Rcpp::Named("sigma")=sigma
	);
}
