// The weighted least squares of LD score regression, shared by the LDSC modules (ldsc_univ.cpp, ldsc_bicov.cpp,
// ldsc_trans.cpp, ldsc_strat.cpp). The SNPs fall into B groups (resampling blocks); a fit only needs the per-group
// cross products X'WX and X'Wy, so every fit (all SNPs, one bootstrap replicate, one leave-one-out set) is a
// weighted sum of the group sums followed by one K x K solve.
//  equal_groups   B contiguous groups of nearly equal size, as LDSC's block bootstrap
//  block_sums     the per-group X'WX (K x K) and X'Wy (K)
//  combine        sum_b c_b (X'WX)_b, sum_b c_b (X'Wy)_b for group weights c (all ones, bootstrap counts, or a
//                 leave-one-out indicator)
//  solve_scaled   A x = b on the unit-diagonal scale, by column-pivoting QR with rank tolerance tol
//  fit_intercept  the solve with the first coefficient (the intercept) estimated freely, fixed at a value, or
//                 estimated inside [lo, hi] (set to the bound it crosses, the others re-solved)
#ifndef PLINK4R_LDSC_KERNEL_H
#define PLINK4R_LDSC_KERNEL_H

#include <RcppEigen.h>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace plink4r {
namespace ldsc {

// 0-based group of each of M rows: nblock groups of floor(M / nblock) rows, the remainder rem either forming one
// more group (rem >= half a group) or added to the last one (LDSC / ldscR rule). Returns the number of groups.
inline int equal_groups(int M, int nblock, std::vector<int>& grp) {
  if (nblock < 1 || M < nblock) throw std::runtime_error("nblock must lie in [1, number of SNPs]");
  const int size = M / nblock, rem = M - size * nblock;
  std::vector<int> len;
  if (rem >= 0.5 * size) {
    len.assign(nblock, size);
    if (rem > 0) len.push_back(rem);
  } else {
    len.assign(nblock, size);
    len.back() += rem;
  }
  grp.resize(M);
  int r = 0;
  for (int b = 0; b < static_cast<int>(len.size()); ++b)
    for (int i = 0; i < len[b]; ++i) grp[r++] = b;
  return static_cast<int>(len.size());
}

struct Sums {
  int K = 0, B = 0;
  std::vector<Eigen::MatrixXd> xtx;
  std::vector<Eigen::VectorXd> xty;
};

// X: M x K (column-major), y, w: M, grp: 0-based group of each row (0..B-1).
inline Sums block_sums(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& w,
                       const std::vector<int>& grp, int B) {
  const int M = static_cast<int>(X.rows()), K = static_cast<int>(X.cols());
  Sums s;
  s.K = K; s.B = B;
  s.xtx.assign(B, Eigen::MatrixXd::Zero(K, K));
  s.xty.assign(B, Eigen::VectorXd::Zero(K));
  for (int i = 0; i < M; ++i) {
    const int b = grp[i];
    const double wi = w[i];
    Eigen::MatrixXd& A = s.xtx[b];
    Eigen::VectorXd& v = s.xty[b];
    for (int c = 0; c < K; ++c) {
      const double xc = wi * X(i, c);
      v[c] += xc * y[i];
      for (int d = 0; d <= c; ++d) A(c, d) += xc * X(i, d);
    }
  }
  for (int b = 0; b < B; ++b) s.xtx[b].triangularView<Eigen::StrictlyUpper>() = s.xtx[b].transpose();
  return s;
}

// sum_b c[b] * group b
inline void combine(const Sums& s, const double* c, Eigen::MatrixXd& A, Eigen::VectorXd& v) {
  A.setZero(s.K, s.K);
  v.setZero(s.K);
  for (int b = 0; b < s.B; ++b) {
    if (c[b] == 0) continue;
    A.noalias() += c[b] * s.xtx[b];
    v.noalias() += c[b] * s.xty[b];
  }
}

inline Eigen::VectorXd solve_scaled(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double tol) {
  const Eigen::VectorXd d = A.diagonal().cwiseSqrt().cwiseInverse();
  const Eigen::MatrixXd As = d.asDiagonal() * A * d.asDiagonal();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(As);
  qr.setThreshold(tol);
  if (qr.rank() < As.cols()) throw std::runtime_error("singular regression: the LD-score columns are collinear");
  return d.cwiseProduct(qr.solve(d.cwiseProduct(b)));
}

// mode 0: intercept estimated freely; 1: fixed at `fixed`; 2: estimated inside [lo, hi].
struct Intercept {
  int mode = 0;
  double fixed = 1, lo = 0, hi = 0;
};

inline Eigen::VectorXd fit_at(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double a, double tol) {
  const int K = static_cast<int>(A.rows());
  Eigen::VectorXd out(K);
  out[0] = a;
  out.tail(K - 1) = solve_scaled(A.bottomRightCorner(K - 1, K - 1), b.tail(K - 1) - a * A.col(0).tail(K - 1), tol);
  return out;
}

inline Eigen::VectorXd fit_intercept(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Intercept& ic,
                                     double tol) {
  if (ic.mode == 1) return fit_at(A, b, ic.fixed, tol);
  Eigen::VectorXd r = solve_scaled(A, b, tol);
  if (ic.mode == 2) {
    if (r[0] < ic.lo) return fit_at(A, b, ic.lo, tol);
    if (r[0] > ic.hi) return fit_at(A, b, ic.hi, tol);
  }
  return r;
}

// One fit per row of counts (R x B, column-major: replicate r uses group b counts(r, b) times), plus a fixed
// penalty P added to every X'WX. Returns R x K.
inline Eigen::MatrixXd refit(const Sums& s, const Rcpp::NumericMatrix& counts, const Eigen::MatrixXd& P,
                             const Intercept& ic, double tol) {
  const int R = counts.nrow();
  if (counts.ncol() != s.B) throw std::runtime_error("the resampling counts do not match the number of blocks");
  Eigen::MatrixXd est(R, s.K);
  std::vector<double> c(s.B);
  Eigen::MatrixXd A;
  Eigen::VectorXd v;
  for (int r = 0; r < R; ++r) {
    for (int b = 0; b < s.B; ++b) c[b] = counts(r, b);
    combine(s, c.data(), A, v);
    est.row(r) = fit_intercept(A + P, v, ic, tol).transpose();
  }
  return est;
}

// Leave-one-group-out fits: B x K.
inline Eigen::MatrixXd leave_one_out(const Sums& s, const Eigen::MatrixXd& P, const Intercept& ic, double tol) {
  std::vector<double> c(s.B, 1.0);
  Eigen::MatrixXd A_all, A_b;
  Eigen::VectorXd v_all;
  combine(s, c.data(), A_all, v_all);
  Eigen::MatrixXd est(s.B, s.K);
  for (int b = 0; b < s.B; ++b)
    est.row(b) = fit_intercept(A_all + P - s.xtx[b], v_all - s.xty[b], ic, tol).transpose();
  return est;
}

// Sample standard deviation of each column (NA with fewer than 2 rows).
inline Rcpp::NumericVector col_sd(const Eigen::MatrixXd& est) {
  const int n = static_cast<int>(est.rows()), K = static_cast<int>(est.cols());
  Rcpp::NumericVector out(K, NA_REAL);
  if (n < 2) return out;
  for (int c = 0; c < K; ++c) out[c] = std::sqrt((est.col(c).array() - est.col(c).mean()).square().sum() / (n - 1));
  return out;
}

// Leave-one-out jackknife standard deviation of each column: sqrt((B - 1) / B * sum (e_b - mean)^2).
inline Rcpp::NumericVector col_jk(const Eigen::MatrixXd& est) {
  const int n = static_cast<int>(est.rows()), K = static_cast<int>(est.cols());
  Rcpp::NumericVector out(K, NA_REAL);
  if (n < 2) return out;
  for (int c = 0; c < K; ++c)
    out[c] = std::sqrt((n - 1.0) / n * (est.col(c).array() - est.col(c).mean()).square().sum());
  return out;
}

// Fit on all rows directly (no grouping).
inline Eigen::VectorXd fit_all(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& w,
                               const Intercept& ic, double tol) {
  const Eigen::MatrixXd Xw = w.asDiagonal() * X;
  const Eigen::MatrixXd A = X.transpose() * Xw;
  const Eigen::VectorXd v = Xw.transpose() * y;
  return fit_intercept(A, v, ic, tol);
}

}  // namespace ldsc
}  // namespace plink4r
#endif
