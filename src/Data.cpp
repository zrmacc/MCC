// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>

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

// ----------------------------------------------------------------------------

// Simulate multivariate recurrent event data for a single subject.
//
// @param censoring_rate Rate of censoring.
// @param death_rate Rate of terminal events.
// @param idx Subject index.
// @param event_rates Vector of recurrent event rates (length K).
// @param tau Truncation time.
// @return Matrix with columns idx, status, time, event_type.

arma::mat SimMvSubjCpp(
	const double censoring_rate,
	const double death_rate,
	const int idx,
	const arma::colvec event_rates,
	const double tau
) {

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

	// Simulate recurrent events for each event type.
	Rcpp::NumericVector id;
	Rcpp::NumericVector status;
	Rcpp::NumericVector time;
	Rcpp::NumericVector event_type;
	const int K = event_rates.n_elem;

	for (int k = 0; k < K; k++) {
		const double event_rate = event_rates(k);
		double follow_up = 0.0;
		while (follow_up <= obs_time) {
			double gap = (-1.0 / event_rate) * std::log(arma::randu());
			follow_up += gap;
			if (follow_up <= obs_time) {
				id.push_back(idx);
				status.push_back(1.0);
				time.push_back(follow_up);
				event_type.push_back(static_cast<double>(k + 1));
			}
		}
	}

	// Terminal row.
	id.push_back(idx);
	status.push_back(final_status);
	time.push_back(obs_time);
	event_type.push_back(Rcpp::NumericVector::get_na());

	// Sort by time, then event_type (NA last within ties).
	const int n_rows = id.size();
	arma::uvec order(n_rows);
	std::iota(order.begin(), order.end(), 0);
	std::stable_sort(
		order.begin(),
		order.end(),
		[&](const arma::uword a, const arma::uword b) {
			if (time[a] < time[b]) {
				return true;
			}
			if (time[a] > time[b]) {
				return false;
			}
			const bool a_na = Rcpp::NumericVector::is_na(event_type[a]);
			const bool b_na = Rcpp::NumericVector::is_na(event_type[b]);
			if (a_na && !b_na) {
				return false;
			}
			if (!a_na && b_na) {
				return true;
			}
			if (a_na && b_na) {
				return false;
			}
			return event_type[a] < event_type[b];
		}
	);

	Rcpp::NumericVector id_out(n_rows);
	Rcpp::NumericVector status_out(n_rows);
	Rcpp::NumericVector time_out(n_rows);
	Rcpp::NumericVector event_type_out(n_rows);
	for (int r = 0; r < n_rows; r++) {
		const arma::uword j = order(r);
		id_out[r] = id[j];
		status_out[r] = status[j];
		time_out[r] = time[j];
		event_type_out[r] = event_type[j];
	}

	arma::mat out = arma::join_horiz(
		arma::colvec(id_out),
		arma::colvec(status_out),
		arma::colvec(time_out),
		arma::colvec(event_type_out)
	);
	return out;
}

// ----------------------------------------------------------------------------

//' Simulate multivariate recurrent event data for multiple subjects.
//'
//' @param censoring_rate Rate of censoring.
//' @param death_rate Rate of terminal events.
//' @param idx Subject index.
//' @param event_rates Matrix of recurrent event rates (n subjects x K types).
//' @param tau Truncation time.
//' @return Long-format recurrent event data with event_type.
// [[Rcpp::export]]

SEXP SimMvDataCpp(
	const arma::colvec censoring_rate,
	const arma::colvec death_rate,
	const arma::colvec idx,
	const arma::mat event_rates,
	const double tau
) {

	arma::mat subj;
	arma::mat out;
	for (int i = 0; i < idx.n_elem; i++) {
		subj = SimMvSubjCpp(
			censoring_rate(i),
			death_rate(i),
			static_cast<int>(idx(i)),
			event_rates.row(i).t(),
			tau
		);
		out = arma::join_vert(out, subj);
	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("idx") = out.col(0),
		Rcpp::Named("status") = out.col(1),
		Rcpp::Named("time") = out.col(2),
		Rcpp::Named("event_type") = out.col(3)
	);
}


