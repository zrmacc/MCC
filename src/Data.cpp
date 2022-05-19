// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------

// Simulate Data for a Single Subject
// 
// @param censoring_rate Rate of censoring. 
// @param death_rate Rate of terminal events. 
// @param idx Subject index.
// @param event_rate Rate of events.
// @param tau Truncation time.
// @return Recurrent event data for a single subject.

arma::mat SimSubjCpp(
	const double censoring_rate,
	const double death_rate,
	const int idx,
	const double event_rate,
	const double tau
){

	// Censoring.
	double cens;
	if (censoring_rate > 0.0) {
		cens = (-1.0 / censoring_rate) * std::log(arma::randu());
	  cens = std::min(cens, tau);
	} else {
		cens = tau;
	}

	// Death.
	double death;
	if (death_rate > 0.0) {
		death = (-1.0 / death_rate) * std::log(arma::randu());
	} else {
		death = std::numeric_limits<double>::infinity();
	}
	
	// Final status.
	double final_status;
	if (death <= cens) {
	  final_status = 2.0;
	} else {
	  final_status = 0.0;
	}
	const double obs_time = std::min(cens, death);
	
	// Simulate event times.
	Rcpp::NumericVector id;
	Rcpp::NumericVector status;
	Rcpp::NumericVector time;
	double follow_up = 0.0;
	while (follow_up <= obs_time) {
	  
	  // Gap time.
	  double gap = (-1.0 / event_rate) * std::log(arma::randu());
	  follow_up += gap;
	  
	  // If follow-up <= observation time, append event.
	  if (follow_up <= obs_time) {
	    id.push_back(idx);
	    status.push_back(1.0);
	    time.push_back(follow_up);
	  }
	  
	}
	
	// Add final time.
	id.push_back(idx);
	status.push_back(final_status);
	time.push_back(obs_time);

	// Convert to matrix.
	arma::mat out = arma::join_horiz(
		arma::colvec(id),
		arma::colvec(status),
		arma::colvec(time)
	);
	return out;
}

// ----------------------------------------------------------------------------

//' Simulate Data for Multiple Subjects
//' 
//' @param censoring_rate Rate of censoring. 
//' @param death_rate Rate of terminal events. 
//' @param idx Subject index.
//' @param event_rate Rate of events.
//' @param tau Truncation time.
//' @return Recurrent event data for a single subject.
// [[Rcpp::export]]

SEXP SimDataCpp(
	const arma::colvec censoring_rate,
	const arma::colvec death_rate,
	const arma::colvec idx,
	const arma::colvec event_rate,
	const double tau
){

	// Loop over subjects.
	arma::mat subj;
	arma::mat out;
	for(int i=0; i<idx.size(); i++) {

		subj = SimSubjCpp(
			censoring_rate(i),
			death_rate(i),
			idx(i),
			event_rate(i),
			tau
		);

		out = arma::join_vert(out, subj);

	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("idx")=out.col(0),
		Rcpp::Named("status")=out.col(1),
		Rcpp::Named("time")=out.col(2)
	);
}


