// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------

//' Calculate AUC Influence Function Contributions 
//'
//' @param event_rate Event rate. 
//' @param idx Unique subject index. 
//' @param haz Hazard.
//' @param nar Number at risk.
//' @param status Status indicator for a single subject.
//' @param surv Survival. 
//' @param time Observation times for a single subject.
//' @return Numeric variance.
//' @noRd 
// [[Rcpp::export]]

SEXP PsiAUC(
    const arma::colvec event_rate,
    const arma::colvec idx,
    const arma::colvec haz,
    const arma::colvec nar,
    const arma::colvec status,
    const arma::colvec surv,
    const double tau,
    const arma::colvec time
){

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  const arma::colvec unique_times = arma::unique(time);
  const int n_times = unique_times.size();

  // Loop over subjects.
  Rcpp::NumericVector id;
  Rcpp::NumericVector influence;
  for(int i=0; i<n; i++) {
    
    // Current subject.
    const int key = unique_idx(i);
    
    // Subject's time and status.
    const arma::colvec subj_time = time.elem(arma::find(idx == key));
    const arma::colvec subj_status = status.elem(arma::find(idx == key));
    
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
      
      // Check if subject is censored.
      if (arma::any(subj_time == current_time && subj_status == 0.0)) {
        at_risk = 0.0;
      }			
      
      // Check if subject has an event.
      if (arma::any(subj_time == current_time && subj_status == 1.0)) {
        event(j) = 1.0;
      }
      
      // Check if subject dies.
      if (arma::any(subj_time == current_time && subj_status == 2.0)) {
        death(j) = 1.0;
        at_risk = 0.0;
      }
      
    }
    
    // Event martingale.
    const arma::colvec dm_event = event - ind_risk % event_rate;
    
    // Death martingale.
    const arma::colvec dm_death = death - ind_risk % haz;

    // Proportion at risk.
    const arma::colvec prop_risk = nar / nar(0);

    // Influence function.
    const arma::colvec nu = arma::sum((tau - unique_times) % surv % event_rate) - \
    	arma::cumsum((tau - unique_times) % surv % event_rate);
    
    const double psi = arma::sum((tau - unique_times) % surv / prop_risk % dm_event) - \
    	arma::sum(nu / prop_risk % dm_death);

    // Cache.
    id.push_back(key);
    influence.push_back(psi);
	}

	// Output.
	return Rcpp::DataFrame::create(
		Rcpp::Named("idx")=id,
		Rcpp::Named("psi")=influence
	);
}