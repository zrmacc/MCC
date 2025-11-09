// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

// Calculate Variance of MCF
// 
// @param idx Unique subject index. 
// @param event_rate Event rate. 
// @param haz Hazard.
// @param mcf Mean cumulative function.
// @param prop_risk Proportion at risk.
// @param status Status indicators.
// @param surv Survival. 
// @param time Observation times.
// @param weights Jump weights.
// @return Numeric vector of estimated variances at each distinct time.
arma::colvec VarMCF(
    const arma::colvec idx,
    const arma::colvec event_rate,
    const arma::colvec haz,
    const arma::colvec mcf,
    const arma::colvec prop_risk,
    const arma::colvec status,
    const arma::colvec surv, 
    const arma::colvec time,
    const arma::colvec weights
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  const arma::colvec unique_times = arma::unique(time);
  const int n_times = unique_times.size();
  
  // Variance vector.
  arma::colvec var = arma::zeros(n_times);

  // Left-continuous versions for integration.
  arma::colvec surv_left = arma::shift(surv, 1); 
  surv_left(0) = 1.0;

  arma::colvec mcf_left  = arma::shift(mcf, 1); 
  mcf_left(0)  = 0.0;

  // Safe version of proportion at risk:
  // Summand is set to zero where proportion at risk is zero.
  arma::colvec prop_risk_safe = prop_risk;
  prop_risk_safe.elem(arma::find(prop_risk_safe == 0.0)).fill(arma::datum::inf);

  // Time -> index map.
  std::unordered_map<double, arma::uword> time_pos;
  time_pos.reserve(static_cast<size_t>(n_times));
  for (arma::uword j = 0; j < static_cast<arma::uword>(n_times); ++j) {
    time_pos.emplace(unique_times(j), j);
  }

  // Initialize.
  arma::colvec event(n_times);
  arma::colvec death(n_times);
  arma::colvec ind_risk(n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {

    // Zero event, death, and ind_risk.
    event.zeros();
    death.zeros();
    ind_risk.zeros();
    
    // Current subject.
    const int key = unique_idx(i);
    const arma::uvec indices = arma::find(idx == key);
    
    // Subject's time and status.
    const arma::colvec subj_time = time.elem(indices);
    const arma::colvec subj_status = status.elem(indices);
    const arma::colvec subj_weight = weights.elem(indices);
    
    // Find index for last time the subject was at risk.
    const double t_last = subj_time(subj_time.n_elem - 1);
    const arma::uword j_star = time_pos[t_last];
    ind_risk.subvec(0, j_star).ones();

    // Find positions of subject's events.
    for (arma::uword r = 0; r < subj_status.n_elem; ++r) {
      if (subj_status(r) == 1.0) {
        const arma::uword k = time_pos[subj_time(r)];
        event(k) += subj_weight(r);
      }
    }

    // Set death, if present.
    for (arma::uword r = 0; r < subj_status.n_elem; ++r) {
      if (subj_status(r) == 2.0) {
        const arma::uword k = time_pos[subj_time(r)];
        death(k) = 1.0;
        break;
      }
    }
    
    // Event martingale.
    const arma::colvec dm_event = event - ind_risk % event_rate;
    
    // Death martingale.
    const arma::colvec dm_death = death - ind_risk % haz;

    // Influence function.
    const arma::colvec psi = arma::cumsum(surv_left / prop_risk_safe % dm_event) \
      - mcf % arma::cumsum(dm_death / prop_risk_safe)                       \
      + arma::cumsum(mcf_left % dm_death / prop_risk_safe);
      
      // Variance contribution.
      var += (psi % psi) / n;
  }

  // Enforce monotonicity.
  arma::colvec var_mon = var;
  double run_max = var_mon(0);
  for (arma::uword j = 1; j < var_mon.n_elem; ++j) {
    if (var_mon(j) < run_max) var_mon(j) = run_max;
    else run_max = var_mon(j);
  }

  return var_mon;
}


// ----------------------------------------------------------------------------

//' Calculate Mean Cumulative Function
//' 
//' Tabulates the mean cumulative function. See equation 2.1 of 
//'  <doi:10.1111/j.0006-341X.2000.00554.x>.
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
//' @param time Observation time.
//' @param weights Jump weights.
//' @param calc_var Calculate variance of the MCF?
//' @return Data.frame with these columns:
//' \itemize{
//'    \item `times`, distinct observation times.
//'    \item `censor`, number of censorings.
//'    \item `death`, number of deaths.
//'    \item `event`, number of events.
//'    \item `haz`, instantaneous hazard (of death).
//'    \item `surv`, survival probability.
//'    \item `event_rate`, instantaneous event rate.
//'    \item `mcf`, mean cumulative function. 
//'    \item `se_mcf`, standard error of the MCF.
//' }
//' @noRd
// [[Rcpp::export]]
SEXP CalcMCFCpp(
	const arma::colvec idx,
	const arma::colvec status,
	const arma::colvec time,
	const arma::colvec weights,
	const bool calc_var = true
){

	// Tabulte the number of censorings, events, and deaths,
	// at each unique event time.

	// Subjects.
	const arma::colvec unique_idx = arma::unique(idx);
	const int n = unique_idx.size();

	// Unique times.
	const arma::colvec unique_times = arma::unique(time);
	const int n_times = unique_times.size();

	// Censoring, death, event, and NAR counts:
	// Length = number of unique times.
	arma::colvec censor(n_times);
	arma::colvec death(n_times);
	arma::colvec event(n_times);
	arma::colvec event_weighted(n_times);
	arma::colvec nar(n_times);

	// Loop over unique times.
	double current_nar = n;
	for(int i=0; i<n_times; i++) {

		// Store initial number at risk.
		nar(i) = current_nar;

		// Logical vector of indices for the current event time.
		arma::uvec indices = arma::find(time == unique_times(i));

    // Get status and weights at the indices.
    arma::vec current_status = status.elem(indices);
    arma::vec current_weights = weights.elem(indices);

		// Count occurrences of each event type.
    int n_zeros = arma::accu(current_status == 0.0);
    int n_ones = arma::accu(current_status == 1.0);
    int n_twos = arma::accu(current_status == 2.0);

    // Store event types.
    censor(i) = n_zeros;
    event(i) = n_ones;
    death(i) = n_twos;

    // Store weighted event sum.
    event_weighted(i) = arma::accu(current_weights % (current_status == 1.0));

		// Update NAR.
		current_nar -= n_zeros + n_twos;
	}

	// Hazard (of death).
	const arma::colvec haz = death / nar;

	// Survival probability.
	const arma::colvec surv = arma::cumprod(1 - haz);

	// Left continuous version of survival function.
	arma::colvec surv_left = arma::shift(surv, 1);
	surv_left(0) = 1.0;

	// Event rate.
	const arma::colvec event_rate = event / nar;

	// Weighted event rate.
	const arma::colvec weighted_event_rate = event_weighted / nar;

	// Mean cumulative function.
	const arma::colvec mcf = arma::cumsum(surv_left % weighted_event_rate);

	// Proportion at risk.
	const arma::colvec prop_risk = nar / nar(0);

	// Variance vector.
	arma::colvec var = arma::zeros(n_times);

	// Calculate variance.
	if (calc_var) {

		var = VarMCF(idx, weighted_event_rate, haz, mcf, 
			prop_risk, status, surv, time, weights);

	}

	// Standard error
	const arma::colvec se = arma::sqrt(var / n);

	// Output.
	return Rcpp::DataFrame::create(
		Rcpp::Named("time")=unique_times,
		Rcpp::Named("censor")=censor,
		Rcpp::Named("death")=death,
		Rcpp::Named("events")=event,
		Rcpp::Named("nar")=nar,
		Rcpp::Named("haz")=haz,
		Rcpp::Named("surv")=surv,
		Rcpp::Named("event_rate")=event_rate,
		Rcpp::Named("weighted_event_rate")=weighted_event_rate,
		Rcpp::Named("mcf")=mcf,
		Rcpp::Named("var_mcf")=var,
		Rcpp::Named("se_mcf")=se
	);

}
