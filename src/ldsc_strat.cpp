// Stratified LD score regression on the squared-effect scale (ldscR::ldsc.strat), adapted from plink4r:
//   E[Z_j^2 / N_j] = a / N_j + sum_c tau_c L_jc ,  Z truncated to [-z_max, z_max].
// Regressors (1 / N_j, L_j1..L_jC), L_j1 the overall LD score. Weights
// w_j = 1 / max(1, wL_j) / (1 / N_j + h2 L_j1 / M)^2 with M the panel SNPs and h2 = sum_c tau_c acov_c, acov the
// annotation sums over the whole panel; started at h2 = 0.1 M / M... (w_init) and updated n_reweight times.
// a is estimated inside [lo, hi] or fixed at 1; optional two-step (a from chi^2 < two_step).
// SEs: leave-one-block-out jackknife over `block` (jackknife = TRUE) or bootstrap with `counts` over `block`.
// Ridge: penalty lambda * u * sum_c v_c tau_c^2, lambda given or chosen by leave-one-block-out CV.
// Apanel: the annotations of all panel SNPs (for acov, the overlap matrix A'A and the variances).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ldsc_fit.h"

namespace {
double median(std::vector<double> v) {
  const std::size_t n = v.size();
  std::sort(v.begin(), v.end());
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List ldsc_strat_cpp(const Eigen::Map<Eigen::VectorXd> z, const Eigen::Map<Eigen::VectorXd> N,
                          const Eigen::Map<Eigen::MatrixXd> L, const Eigen::Map<Eigen::VectorXd> wL,
                          const Eigen::Map<Eigen::MatrixXd> Apanel, const Rcpp::IntegerVector& block,
                          bool jackknife, const Rcpp::NumericMatrix& counts, bool intercept, double lo, double hi,
                          double two_step, double z_max, int n_reweight, bool ridge, double lambda,
                          const Rcpp::NumericVector& lambda_grid, double tol) {
  using namespace ldscr;
  const int n = static_cast<int>(z.size()), C = static_cast<int>(L.cols()), K = C + 1;
  if (N.size() != n || L.rows() != n || wL.size() != n || block.size() != n) Rcpp::stop("inputs differ in length");
  if (Apanel.cols() != C) Rcpp::stop("L and the annotations differ in their number of columns");
  if (n <= K) Rcpp::stop("fewer SNPs than regressors");
  std::vector<int> grp(block.begin(), block.end());
  const int B = *std::max_element(grp.begin(), grp.end()) + 1;
  const double M = static_cast<double>(Apanel.rows());

  const Eigen::VectorXd acov = Apanel.colwise().sum().transpose();
  const Eigen::VectorXd amean = acov / M;
  const Eigen::VectorXd avar = ((Apanel.array().square().colwise().sum().transpose() / M) - amean.array().square()).max(0.0).matrix();
  const Eigen::MatrixXd overlap = Apanel.transpose() * Apanel;

  Eigen::MatrixXd X(n, K);
  Eigen::VectorXd y(n), lw(n), invN(n);
  std::vector<char> s1(n);
  int n_trunc = 0;
  for (int j = 0; j < n; ++j) {
    const double az = std::abs(z[j]);
    if (az > z_max) ++n_trunc;
    const double zt = std::min(az, z_max);
    invN[j] = 1.0 / N[j];
    y[j] = zt * zt * invN[j];
    lw[j] = 1.0 / std::max(1.0, wL[j]);
    X(j, 0) = invN[j];
    s1[j] = zt * zt < two_step;
  }
  X.rightCols(C) = L;
  Intercept ic;
  ic.mode = intercept ? 2 : 1;
  ic.fixed = 1.0; ic.lo = lo; ic.hi = hi;
  const bool two = intercept && two_step > 0;

  try {
    auto weights = [&](double h2) {
      return Eigen::VectorXd(lw.array() / (invN.array() + std::max(h2, 0.0) * L.col(0).array() / M).square());
    };
    Eigen::VectorXd w = weights(0.1);
    Eigen::VectorXd b = fit2_rows(X, y, w, s1, ic, two, tol);
    double h2 = b.tail(C).dot(acov);
    for (int it = 0; it < n_reweight; ++it) {
      w = weights(h2);
      b = fit2_rows(X, y, w, s1, ic, two, tol);
      h2 = b.tail(C).dot(acov);
    }

    const Sums2 s = block_sums2(X, y, w, s1, grp, B, two);
    Eigen::MatrixXd Pen = Eigen::MatrixXd::Zero(K, K);
    double lambda_used = NA_REAL;
    Rcpp::NumericVector cv;
    if (ridge) {
      std::vector<double> ones(B, 1.0);
      Eigen::MatrixXd A2, A1;
      Eigen::VectorXd v2, v1;
      combine(s.s2, ones.data(), A2, v2);
      if (two) combine(s.s1, ones.data(), A1, v1);
      Eigen::VectorXd pen_d(K);
      pen_d[0] = 0;
      pen_d.tail(C) = avar;
      std::vector<double> ratio;
      for (int c = 0; c < K; ++c) if (pen_d[c] > 0) ratio.push_back(A2(c, c) / pen_d[c]);
      if (ratio.empty()) Rcpp::stop("ridge: no annotation column varies over the SNPs");
      const double unit = median(ratio);
      if (std::isnan(lambda)) {
        const int G = lambda_grid.size();
        cv = Rcpp::NumericVector(G);
        for (int g = 0; g < G; ++g) {
          const Eigen::MatrixXd Pg = (lambda_grid[g] * unit * pen_d).asDiagonal();
          double e = 0;
          for (int bb = 0; bb < B; ++bb) {
            const Eigen::MatrixXd B1 = two ? Eigen::MatrixXd(A1 + Pg - s.s1.xtx[bb]) : A2;
            const Eigen::VectorXd u1 = two ? Eigen::VectorXd(v1 - s.s1.xty[bb]) : v2;
            const Eigen::VectorXd f = fit2(B1, u1, A2 + Pg - s.s2.xtx[bb], v2 - s.s2.xty[bb], ic, two, tol);
            e += -2.0 * f.dot(s.s2.xty[bb]) + f.dot(s.s2.xtx[bb] * f);
          }
          cv[g] = e;
        }
        const int best = static_cast<int>(std::min_element(cv.begin(), cv.end()) - cv.begin());
        lambda_used = lambda_grid[best];
        cv = cv - cv[best];
      } else {
        lambda_used = lambda;
      }
      Pen = (lambda_used * unit * pen_d).asDiagonal();
      b = fit2(two ? Eigen::MatrixXd(A1 + Pen) : A2, v1, A2 + Pen, v2, ic, two, tol);
      h2 = b.tail(C).dot(acov);
    }

    Eigen::MatrixXd est(0, K);
    if (jackknife) est = leave_one_out2(s, Pen, ic, two, tol);
    else if (counts.nrow() > 0) est = refit2(s, counts, Pen, ic, two, tol);

    return Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(b), Rcpp::Named("h2") = h2, Rcpp::Named("est") = Rcpp::wrap(est),
                              Rcpp::Named("acov") = Rcpp::wrap(acov), Rcpp::Named("avar") = Rcpp::wrap(avar),
                              Rcpp::Named("overlap") = Rcpp::wrap(overlap), Rcpp::Named("lambda") = lambda_used,
                              Rcpp::Named("cv") = cv, Rcpp::Named("n_truncated") = n_trunc, Rcpp::Named("B") = B,
                              Rcpp::Named("M") = M);
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
