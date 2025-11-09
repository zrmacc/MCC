// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Calculate AUC Influence Function Contributions 
//'
//' @param event_rate Event rate. 
//' @param grid_time Times at which the (event_rate, haz, nar, surv) are
//'   evaluated.
//' @param idx Unique subject index. 
//' @param haz Hazard.
//' @param nar Number at risk.
//' @param status Status indicator for a single subject.
//' @param surv Survival. 
//' @param time Observation times for a single subject.
//' @param weights Jump weights.
//' @return Numeric variance.
//' @noRd 
// [[Rcpp::export]]

SEXP PsiAUC(
    const arma::colvec event_rate,
    const arma::colvec grid_time,
    const arma::colvec idx,
    const arma::colvec haz,
    const arma::colvec nar,
    const arma::colvec status,
    const arma::colvec surv,
    const double tau,
    const arma::colvec time,
    const arma::colvec weights
){

  // Unique times.
  const arma::colvec unique_times = arma::unique(time);

  // Confirm that all unique values of time are present in grid_time.
  if (unique_times.n_elem != grid_time.n_elem || !arma::all(unique_times == grid_time)) {
    Rcpp::stop("Unique times do not match grid times.");
  }

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();

  // Left continuous survival function.
  arma::colvec surv_left = arma::shift(surv, 1);
  surv_left(0) = 1.0;

  // Proportion at risk.
  const arma::colvec prop_risk = nar / nar(0);

  // Safe version of proportion at risk:
  // Summand is set to zero where proportion at risk is zero.
  arma::colvec prop_risk_safe = prop_risk;
  prop_risk_safe.elem(arma::find(prop_risk_safe == 0.0)).fill(arma::datum::inf);

  // Clamp temporal difference at zero.
  arma::colvec w = tau - grid_time;
  w = arma::clamp(w, 0.0, arma::datum::inf);

  // nu_k = sum_{m} w_m surv_left_m r_m  -  cumsum_k(w surv_left r)
  arma::colvec nu = arma::sum(w % surv_left % event_rate) - \
    arma::cumsum(w % surv_left % event_rate);

  // a_k = w_k * surv_left_k / prop_risk_k
  arma::colvec a = (w % surv_left) / prop_risk_safe;

  // b_k = nu_k / prop_risk_k
  arma::colvec b = nu / prop_risk_safe;

  // Prefix sums:
  // pref_a[k] = sum_{m<=k} a_m * rate_m.
  arma::colvec pref_a = arma::cumsum(a % event_rate);
  // pref_b[k] = sum_{m<=k} b_m * haz_m.
  arma::colvec pref_b = arma::cumsum(b % haz);

  // Time -> index map.
  const int n_times = grid_time.size();
  std::unordered_map<double, arma::uword> time_pos;
  time_pos.reserve(n_times);
  for (arma::uword j = 0; j < n_times; ++j) {
    time_pos.emplace(grid_time(j), j);
  }

  // Allocate output structures. 
  arma::colvec id = arma::zeros(n);
  arma::colvec influence = arma::zeros(n);

  // Loop over subjects.
  for (int i = 0; i < n; ++i) {

    // Current subject.
    const int key = unique_idx(i);
    const arma::uvec indices = arma::find(idx == key);

    // Subject data.
    const arma::colvec subj_time = time.elem(indices);
    const arma::colvec subj_status = status.elem(indices);
    const arma::colvec subj_weight = weights.elem(indices);

    // Find the index of the subjects last observation time.
    double t_last = subj_time(subj_time.n_elem - 1);
    arma::uword j_star = time_pos[t_last];

    // Sum contributions at the subject's event rows: sum a[k] * weight_k.
    double sum_events = 0.0;
    for (arma::uword r = 0; r < subj_status.n_elem; ++r) {
      if (subj_status(r) == 1.0) {
        arma::uword k = time_pos[subj_time(r)];
        sum_events += a(k) * subj_weight(r);
      }
    }

    // Death jump term: - b[k_death] if the subject dies; 0 otherwise.
    double death_jump = 0.0;
    for (arma::uword r = 0; r < subj_status.n_elem; ++r) {
      if (subj_status(r) == 2.0) {
        arma::uword k = time_pos[subj_time(r)];
        death_jump = b(k);
        break;
      }
    }

    // Influence.
    double psi_i = sum_events - pref_a(j_star) - death_jump + pref_b(j_star);

    // Store.
    id[i] = key;
    influence[i] = psi_i;
  }

	// Output.
	return Rcpp::DataFrame::create(
		Rcpp::Named("idx")=id,
		Rcpp::Named("psi")=influence
	);
}
