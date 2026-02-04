// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------

//' Calculate Augmentation Components
//'
//' @param covars One row per subject.
//' @param psi Influence functions.
//' @return List.
//' @noRd 
// [[Rcpp::export]]
SEXP CalcAugComp(
	const arma::mat covars,
	const arma::colvec psi
){
	// Subjects.
	const int n = covars.n_rows;

	// Residual matrix.
	const arma::rowvec xbar = arma::mean(covars, 0);
	const arma::mat resid = covars.each_row() - xbar;

	// Gamma.
	const arma::colvec gamma = resid.t() * psi / (n * n);

	// Sigma.
	const arma::mat sigma = resid.t() * resid / (n * n);

  // Output.
	return Rcpp::List::create(
		Rcpp::Named("gamma")=gamma,
		Rcpp::Named("n")=n,
		Rcpp::Named("resid")=resid,
		Rcpp::Named("sigma")=sigma,
		Rcpp::Named("xbar")=xbar.t()
	);
}
