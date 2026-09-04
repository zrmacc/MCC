// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::export]]
arma::mat H12MCFCpp(
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const arma::colvec weights,
    const arma::mat design,
    const arma::colvec grid_time,
    const arma::colvec nar,
    const arma::colvec death,
    const arma::colvec event_weighted,
    const arma::colvec surv,
    const double tau,
    const bool integrate
) {
  const arma::colvec ids = arma::unique(idx);
  const arma::uword n = ids.n_elem;
  const arma::uword m = grid_time.n_elem;
  const arma::uword p = design.n_cols;
  if (design.n_rows != n) Rcpp::stop("The design matrix is not aligned with subjects.");

  arma::uword m_eff = 0;
  while (m_eff < m && grid_time(m_eff) <= tau) ++m_eff;
  if (m_eff == 0) return arma::mat(n, p, arma::fill::zeros);

  arma::mat risk(m_eff, n, arma::fill::zeros);
  arma::mat events(m_eff, n, arma::fill::zeros);
  arma::mat deaths(m_eff, n, arma::fill::zeros);
  std::unordered_map<double, arma::uword> time_pos;
  for (arma::uword u = 0; u < m_eff; ++u) time_pos.emplace(grid_time(u), u);

  for (arma::uword i = 0; i < n; ++i) {
    const arma::uvec rows = arma::find(idx == ids(i));
    const double last = time(rows(rows.n_elem - 1));
    for (arma::uword u = 0; u < m_eff; ++u) if (grid_time(u) <= last) risk(u, i) = 1.0;
    for (arma::uword r = 0; r < rows.n_elem; ++r) {
      const arma::uword row = rows(r);
      const auto pos = time_pos.find(time(row));
      if (pos == time_pos.end()) continue;
      if (status(row) == 1.0) events(pos->second, i) += weights(row);
      if (status(row) == 2.0) deaths(pos->second, i) = 1.0;
    }
  }

  const arma::vec y = nar.head(m_eff) / static_cast<double>(n);
  const arma::vec da = event_weighted.head(m_eff) / static_cast<double>(n);
  const arma::vec da_d = death.head(m_eff) / static_cast<double>(n);
  const arma::vec d_l = da_d / y;
  const arma::vec denom = 1.0 - d_l;
  if (arma::any(y <= 0.0)) Rcpp::stop("The empirical at-risk process is zero at or before tau.");
  if (arma::any(denom <= 0.0)) Rcpp::stop("The death risk set is exhausted at or before tau; use an earlier tau.");

  arma::vec s_left(m_eff, arma::fill::ones);
  if (m_eff > 1) s_left.tail(m_eff - 1) = surv.head(m_eff - 1);
  const arma::vec b = s_left / y;
  arma::vec time_weight(m_eff, arma::fill::ones);
  if (integrate) time_weight = arma::clamp(tau - grid_time.head(m_eff), 0.0, arma::datum::inf);

  const arma::mat centered_risk = risk.each_col() - y;
  const arma::mat centered_events = events.each_col() - da;
  const arma::mat centered_deaths = deaths.each_col() - da_d;
  const arma::rowvec mean_w = arma::mean(design, 0);
  const arma::mat agg_risk = risk * design / static_cast<double>(n) - y * mean_w;
  const arma::mat agg_events = events * design / static_cast<double>(n) - da * mean_w;
  const arma::mat agg_deaths = deaths * design / static_cast<double>(n) - da_d * mean_w;

  auto first_terms = [&](const arma::vec& hy, const arma::vec& hd) -> std::pair<arma::vec, arma::vec> {
    const arma::vec dl_h = hd / y - hy % da_d / arma::square(y);
    arma::vec j_left(m_eff, arma::fill::zeros);
    double running = 0.0;
    for (arma::uword u = 0; u < m_eff; ++u) {
      j_left(u) = running;
      running += dl_h(u) / denom(u);
    }
    const arma::vec r_h = j_left + hy / y;
    return std::make_pair(dl_h, r_h);
  };

  arma::mat out(n, p, arma::fill::zeros);
  for (arma::uword k = 0; k < p; ++k) {
    const arma::vec hy = agg_risk.col(k);
    const arma::vec ha = agg_events.col(k);
    const arma::vec hd = agg_deaths.col(k);
    const auto h_terms = first_terms(hy, hd);
    const arma::vec dl_h = h_terms.first;
    const arma::vec r_h = h_terms.second;

    for (arma::uword i = 0; i < n; ++i) {
      const arma::vec ky = centered_risk.col(i);
      const arma::vec ka = centered_events.col(i);
      const arma::vec kd = centered_deaths.col(i);
      const auto k_terms = first_terms(ky, kd);
      const arma::vec dl_k = k_terms.first;
      const arma::vec r_k = k_terms.second;
      const arma::vec dl_hk = -hy % kd / arma::square(y) - ky % hd / arma::square(y)
        + 2.0 * hy % ky % da_d / arma::pow(y, 3.0);
      arma::vec j_hk_left(m_eff, arma::fill::zeros);
      double running = 0.0;
      for (arma::uword u = 0; u < m_eff; ++u) {
        j_hk_left(u) = running;
        running += dl_hk(u) / denom(u) + dl_h(u) * dl_k(u) / std::pow(denom(u), 2.0);
      }
      const arma::vec r_hk = j_hk_left - hy % ky / arma::square(y);
      const arma::vec increment = b % (r_h % r_k - r_hk) % da
        - b % r_h % ka - b % r_k % ha;
      out(i, k) = arma::sum(time_weight % increment);
    }
  }
  return out;
}
