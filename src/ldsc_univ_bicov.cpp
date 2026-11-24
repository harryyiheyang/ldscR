// Univariate and bivariate LD score regression (ldscR::ldsc.univ, ldscR::ldsc.bicov).
// univ:  y_j = min(|z_j|, z_max)^2 on (1, l_j), l_j = L_j N_j / M, M = SNPs of the LD panel;
//        w_j = 1 / max(1, wL_j) / (1 + h2 l_j)^2, h2 started at 0.1 and updated n_reweight times.
// bicov: y_j = z1_j z2_j (both truncated) on (1, l_j), l_j = L_j sqrt(N1_j N2_j) / M;
//        w_j = 1 / max(1, wL_j) / (v1_j v2_j + mu_j^2), v_i = 1 + L_j N_ij h2_i / M, mu_j = ecov + gcov l_j
//        (Var(z1 z2) = v1 v2 + mu^2 for bivariate normal z), updated n_reweight times.
// Two-step (two_step > 0): intercept from SNPs with chi^2 < two_step (both traits for bicov), slope from all SNPs.
// Standard errors: block bootstrap over contiguous SNP groups with the final weights; replicate r uses group b
// counts(r, b) times. bicov fits both univariate regressions on the same SNPs and the same draws, so rg has a
// bootstrap SE without a separate function.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include "ldsc_fit.h"

namespace {
using namespace ldscr;

struct UnivOut {
  Eigen::VectorXd b;    // intercept, h2
  Eigen::MatrixXd est;  // R x 2
};

UnivOut univ_core(const Eigen::VectorXd& z, const Eigen::VectorXd& N, const Eigen::VectorXd& L,
                  const Eigen::VectorXd& lw, double M, double z_max, double two_step, int n_reweight,
                  const std::vector<int>& grp, int B, const Rcpp::NumericMatrix& counts, double tol) {
  const int n = static_cast<int>(z.size());
  const bool two = two_step > 0;
  Eigen::MatrixXd X(n, 2);
  Eigen::VectorXd y(n);
  std::vector<char> s1(n);
  for (int j = 0; j < n; ++j) {
    const double zt = std::min(std::abs(z[j]), z_max);
    y[j] = zt * zt;
    X(j, 0) = 1.0;
    X(j, 1) = L[j] * N[j] / M;
    s1[j] = y[j] < two_step;
  }
  const Intercept ic;
  Eigen::VectorXd w = lw.array() / (1.0 + 0.1 * X.col(1).array()).square();
  Eigen::VectorXd b = fit2_rows(X, y, w, s1, ic, two, tol);
  for (int it = 0; it < n_reweight; ++it) {
    w = lw.array() / (1.0 + std::max(b[1], 0.0) * X.col(1).array()).square();
    b = fit2_rows(X, y, w, s1, ic, two, tol);
  }
  UnivOut out;
  out.b = b;
  out.est = Eigen::MatrixXd(0, 2);
  if (counts.nrow() > 0) out.est = refit2(block_sums2(X, y, w, s1, grp, B, two), counts, Eigen::MatrixXd::Zero(2, 2), ic, two, tol);
  return out;
}

Eigen::VectorXd over_count(const Eigen::VectorXd& wL) { return (1.0 / wL.array().max(1.0)).matrix(); }
}  // namespace

// z, N, L, wL: regression SNPs in genomic order; M: panel SNPs; counts: R x B (R = 0: no bootstrap).
// [[Rcpp::export]]
Rcpp::List ldsc_univ_cpp(const Eigen::Map<Eigen::VectorXd> z, const Eigen::Map<Eigen::VectorXd> N,
                         const Eigen::Map<Eigen::VectorXd> L, const Eigen::Map<Eigen::VectorXd> wL, double M,
                         double z_max, double two_step, int nblock, const Rcpp::NumericMatrix& counts,
                         int n_reweight, double tol) {
  const int n = static_cast<int>(z.size());
  if (N.size() != n || L.size() != n || wL.size() != n) Rcpp::stop("inputs differ in length");
  if (n < 3) Rcpp::stop("fewer than 3 SNPs");
  try {
    std::vector<int> grp;
    const int B = plink4r::ldsc::equal_groups(n, std::min(nblock, n), grp);
    const UnivOut u = univ_core(z, N, L, over_count(wL), M, z_max, two_step, n_reweight, grp, B, counts, tol);
    return Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(u.b), Rcpp::Named("boot") = Rcpp::wrap(u.est),
                              Rcpp::Named("B") = B);
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}

// h21, h22 are always estimated here (univariate LDSC on the same SNPs and the same bootstrap draws).
// [[Rcpp::export]]
Rcpp::List ldsc_bicov_cpp(const Eigen::Map<Eigen::VectorXd> z1, const Eigen::Map<Eigen::VectorXd> z2,
                          const Eigen::Map<Eigen::VectorXd> N1, const Eigen::Map<Eigen::VectorXd> N2,
                          const Eigen::Map<Eigen::VectorXd> L, const Eigen::Map<Eigen::VectorXd> wL, double M,
                          double z_max, double two_step, int nblock,
                          const Rcpp::NumericMatrix& counts, int n_reweight, double tol) {
  const int n = static_cast<int>(z1.size());
  if (z2.size() != n || N1.size() != n || N2.size() != n || L.size() != n || wL.size() != n)
    Rcpp::stop("inputs differ in length");
  if (n < 3) Rcpp::stop("fewer than 3 SNPs");
  try {
    std::vector<int> grp;
    const int B = plink4r::ldsc::equal_groups(n, std::min(nblock, n), grp);
    const Eigen::VectorXd lw = over_count(wL);
    const UnivOut u1 = univ_core(z1, N1, L, lw, M, z_max, two_step, 3, grp, B, counts, tol);  // as ldsc.univ
    const UnivOut u2 = univ_core(z2, N2, L, lw, M, z_max, two_step, 3, grp, B, counts, tol);
    const double h21 = u1.b[1], h22 = u2.b[1];
    const bool two = two_step > 0;
    Eigen::MatrixXd X(n, 2);
    Eigen::VectorXd y(n), v(n);
    std::vector<char> s1(n);
    for (int j = 0; j < n; ++j) {
      const double a1 = std::min(std::abs(z1[j]), z_max), a2 = std::min(std::abs(z2[j]), z_max);
      y[j] = (z1[j] < 0 ? -a1 : a1) * (z2[j] < 0 ? -a2 : a2);
      X(j, 0) = 1.0;
      X(j, 1) = L[j] * std::sqrt(N1[j]) * std::sqrt(N2[j]) / M;
      v[j] = (1.0 + L[j] * N1[j] / M * std::max(h21, 0.0)) * (1.0 + L[j] * N2[j] / M * std::max(h22, 0.0));
      s1[j] = a1 * a1 < two_step && a2 * a2 < two_step;
    }
    const Intercept ic;
    Eigen::VectorXd w = lw.cwiseQuotient(v);
    Eigen::VectorXd b = fit2_rows(X, y, w, s1, ic, two, tol);
    for (int it = 0; it < n_reweight; ++it) {
      const Eigen::ArrayXd mu = b[0] + b[1] * X.col(1).array();
      w = lw.array() / (v.array() + mu.square());
      b = fit2_rows(X, y, w, s1, ic, two, tol);
    }
    Eigen::MatrixXd est(0, 2);
    if (counts.nrow() > 0) est = refit2(block_sums2(X, y, w, s1, grp, B, two), counts, Eigen::MatrixXd::Zero(2, 2), ic, two, tol);
    return Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(b), Rcpp::Named("boot") = Rcpp::wrap(est),
                              Rcpp::Named("h21") = h21, Rcpp::Named("h22") = h22,
                              Rcpp::Named("boot1") = Rcpp::wrap(u1.est), Rcpp::Named("boot2") = Rcpp::wrap(u2.est),
                              Rcpp::Named("univ1") = Rcpp::wrap(u1.b), Rcpp::Named("univ2") = Rcpp::wrap(u2.b),
                              Rcpp::Named("B") = B);
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}

// Trans-ancestry bivariate LDSC (ldscR::ldsc.trans): as bicov with the cross-ancestry LD score L12 in the
// regressor and the ancestry-specific scores L1, L2 in v; no over-counting weight (as the original).
// [[Rcpp::export]]
Rcpp::List ldsc_trans_cpp(const Eigen::Map<Eigen::VectorXd> z1, const Eigen::Map<Eigen::VectorXd> z2,
                          const Eigen::Map<Eigen::VectorXd> N1, const Eigen::Map<Eigen::VectorXd> N2,
                          const Eigen::Map<Eigen::VectorXd> L1, const Eigen::Map<Eigen::VectorXd> L2,
                          const Eigen::Map<Eigen::VectorXd> L12, double M, double h21, double h22, int nblock,
                          const Rcpp::NumericMatrix& counts, int n_reweight, double tol) {
  const int n = static_cast<int>(z1.size());
  try {
    std::vector<int> grp;
    const int B = plink4r::ldsc::equal_groups(n, std::min(nblock, n), grp);
    Eigen::MatrixXd X(n, 2);
    Eigen::VectorXd y(n), v(n);
    std::vector<char> s1(n, 1);
    for (int j = 0; j < n; ++j) {
      y[j] = z1[j] * z2[j];
      X(j, 0) = 1.0;
      X(j, 1) = L12[j] * std::sqrt(N1[j]) * std::sqrt(N2[j]) / M;
      v[j] = (1.0 + L1[j] * N1[j] * h21 / M) * (1.0 + L2[j] * N2[j] * h22 / M);
    }
    const Intercept ic;
    Eigen::VectorXd w = v.cwiseInverse();
    Eigen::VectorXd b = fit2_rows(X, y, w, s1, ic, false, tol);
    for (int it = 0; it < n_reweight; ++it) {
      const Eigen::ArrayXd mu = b[0] + b[1] * X.col(1).array();
      w = 1.0 / (v.array() + mu.square());
      b = fit2_rows(X, y, w, s1, ic, false, tol);
    }
    Eigen::MatrixXd est(0, 2);
    if (counts.nrow() > 0) est = refit2(block_sums2(X, y, w, s1, grp, B, false), counts, Eigen::MatrixXd::Zero(2, 2), ic, false, tol);
    return Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(b), Rcpp::Named("boot") = Rcpp::wrap(est), Rcpp::Named("B") = B);
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
