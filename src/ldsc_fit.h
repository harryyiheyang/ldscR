// Two-step fits and resampling on top of ldsc_kernel.h (from plink4r).
// A fit uses two sets of group sums over all SNPs: s1 with the step-1 outcome y1 (the statistic truncated at the
// two-step cutoff; no SNP is dropped) and s2 with the outcome y. With the two-step, the intercept (first
// coefficient) comes from s1 and the slopes from s2 with the intercept fixed; without it, s2 alone.
#ifndef LDSCR_LDSC_FIT_H
#define LDSCR_LDSC_FIT_H

#include "ldsc_kernel.h"

namespace ldscr {
using namespace plink4r::ldsc;

inline Eigen::VectorXd fit2(const Eigen::MatrixXd& A1, const Eigen::VectorXd& v1, const Eigen::MatrixXd& A2,
                            const Eigen::VectorXd& v2, const Intercept& ic, bool two, double tol) {
  if (!two || ic.mode == 1) return fit_intercept(A2, v2, ic, tol);
  const double a = fit_intercept(A1, v1, ic, tol)[0];
  return fit_at(A2, v2, a, tol);
}

inline Eigen::VectorXd fit2_rows(const Eigen::MatrixXd& X, const Eigen::VectorXd& y1, const Eigen::VectorXd& y,
                                 const Eigen::VectorXd& w, const Intercept& ic, bool two, double tol) {
  const Eigen::MatrixXd Xw = w.asDiagonal() * X;
  const Eigen::MatrixXd A2 = X.transpose() * Xw;
  const Eigen::VectorXd v2 = Xw.transpose() * y;
  if (!two) return fit_intercept(A2, v2, ic, tol);
  return fit2(A2, Xw.transpose() * y1, A2, v2, ic, two, tol);
}

struct Sums2 {
  Sums s1, s2;
};

inline Sums2 block_sums2(const Eigen::MatrixXd& X, const Eigen::VectorXd& y1, const Eigen::VectorXd& y,
                         const Eigen::VectorXd& w, const std::vector<int>& grp, int B, bool two) {
  Sums2 out;
  out.s2 = block_sums(X, y, w, grp, B);
  if (two) out.s1 = block_sums(X, y1, w, grp, B);
  return out;
}

// One fit per row of counts (R x B); P is added to both X'WX.
inline Eigen::MatrixXd refit2(const Sums2& s, const Rcpp::NumericMatrix& counts, const Eigen::MatrixXd& P,
                              const Intercept& ic, bool two, double tol) {
  const int R = counts.nrow(), B = s.s2.B, K = s.s2.K;
  if (counts.ncol() != B) throw std::runtime_error("the resampling counts do not match the number of blocks");
  Eigen::MatrixXd est(R, K), A1, A2;
  Eigen::VectorXd v1, v2;
  std::vector<double> c(B);
  for (int r = 0; r < R; ++r) {
    for (int b = 0; b < B; ++b) c[b] = counts(r, b);
    combine(s.s2, c.data(), A2, v2);
    if (two) combine(s.s1, c.data(), A1, v1);
    est.row(r) = fit2(two ? Eigen::MatrixXd(A1 + P) : A2, v1, A2 + P, v2, ic, two, tol).transpose();
  }
  return est;
}

// Leave-one-group-out fits (B x K).
inline Eigen::MatrixXd leave_one_out2(const Sums2& s, const Eigen::MatrixXd& P, const Intercept& ic, bool two,
                                      double tol) {
  const int B = s.s2.B, K = s.s2.K;
  std::vector<double> ones(B, 1.0);
  Eigen::MatrixXd A1, A2;
  Eigen::VectorXd v1, v2;
  combine(s.s2, ones.data(), A2, v2);
  if (two) combine(s.s1, ones.data(), A1, v1);
  Eigen::MatrixXd est(B, K);
  for (int b = 0; b < B; ++b) {
    Eigen::MatrixXd B1 = A2, B2 = A2 + P - s.s2.xtx[b];
    Eigen::VectorXd u1 = v2, u2 = v2 - s.s2.xty[b];
    if (two) { B1 = A1 + P - s.s1.xtx[b]; u1 = v1 - s.s1.xty[b]; }
    est.row(b) = fit2(B1, u1, B2, u2, ic, two, tol).transpose();
  }
  return est;
}

}  // namespace ldscr
#endif
