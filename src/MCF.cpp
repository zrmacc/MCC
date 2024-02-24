// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
  
  // Calculation of martingales. 
  // Note: ensure consistency with PsiAUC.cpp
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    const int key = unique_idx(i);
    
    // Subject's time and status.
    const arma::uvec indices = arma::find(idx == key);
    const arma::colvec subj_time = time.elem(indices);
    const arma::colvec subj_status = status.elem(indices);
    const arma::colvec subj_weight = weights.elem(indices);

    const arma::uvec subj_status_0 = (subj_status == 0.0);
    const arma::uvec subj_status_1 = (subj_status == 1.0);
    const arma::uvec subj_status_2 = (subj_status == 2.0);
    
    // Subject-specific event and status indicators.
    arma::colvec event = arma::zeros(n_times);
    arma::colvec death = arma::zeros(n_times);
    
    // At risk indicator.
    arma::colvec ind_risk = arma::zeros(n_times);
    double at_risk = 1.0;
    
    // Loop over unique times.
    for(int j=0; j<n_times; j++) {
      
      // Current time.
      const double current_time = unique_times(j);
      
      // Risk status.
      ind_risk(j) = at_risk;
      if (at_risk == 0.0) {
        break;
      }

      // Subject time equatls current time.
      const arma::uvec obs_time = (subj_time == current_time);
      
      // Check if subject is censored.
      if (arma::any(obs_time && subj_status_0)) {
        at_risk = 0.0;
      }     
      
      // Check if subject has an event.
      const arma::uvec event_time = (obs_time && subj_status_1);
      if (arma::any(event_time)) {
        event(j) = arma::accu(subj_weight % event_time);
      }
      
      // Check if subject dies.
      if (arma::any(obs_time && subj_status_2)) {
        death(j) = 1.0;
        at_risk = 0.0;
      }
      
    }
    
    // Event martingale.
    const arma::colvec dm_event = event - ind_risk % event_rate;
    
    // Death martingale.
    const arma::colvec dm_death = death - ind_risk % haz;
    
    // Influence function.
    const arma::colvec psi = arma::cumsum(surv / prop_risk % dm_event) \
      - mcf % arma::cumsum(dm_death / prop_risk)                       \
      + arma::cumsum(mcf % dm_death / prop_risk);
      
      // Variance contribution.
      var += (psi % psi) / n;
  }
  return var;
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

	// Event rate.
	const arma::colvec event_rate = event / nar;

	// Weighted event rate.
	const arma::colvec weighted_event_rate = event_weighted / nar;

	// Mean cumulative function.
	const arma::colvec mcf = arma::cumsum(surv % weighted_event_rate);

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
