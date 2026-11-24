// LDAK / SumHer in C++.
//
// ldak_tagging_cpp: for each eigen LD block (eigen_bin.h, from plink4r; SBayesRC or plink4r layout), forms the
// correlations of `chunk` rows at a time, r_jk = u_j' diag(lam) u_k (divided by sqrt(d_j d_k), d = diag of the
// stored factor, when normalize = TRUE), and
//   thin (weight_mode 1): scanning j in order, j gets weight 0 if a later SNP k still kept has r_jk^2 >= thin_r2
//                          (LDAK --thin --window-prune), else 1;
//   q_kc = A_kc w_k [2 f_k (1 - f_k)]^(1 + alpha)                (LDAK expected heritability)
//   v_jc = q_jc + sum_{k != j} r~_jk^2 q_kc,  t_j = reg_j + sum_{k != j} r~_jk^2 reg_k,
// with r~^2 = r^2 - (1 - r^2) / (n_ref - 2) (no correction when n_ref <= 2). The direct chunk product is used
// (single-precision products through the run-time SIMD kernel, la_dispatch.h)
// rather than plink4r's G = sum_k a_k u_k u_k' trick: with C categories it costs m^2 (k + C) per block instead of
// C m k^2, and it gives the r^2 needed by the thinning anyway. OpenMP over blocks, largest first.
//
// sumher_cpp: SumHer's likelihood (optionally ridge-penalised, see sumher_cpp), S_j ~ E_j chi^2_1 with E_j = 1 + T_j' theta, T_j = (n_j / s) v_j (and n_j / s
// for the intercept), each SNP weighted 1 / t_j:
//   l = -1/2 sum_j (S_j / E_j + log(S_j E_j)) / t_j
// maximised by Newton-Raphson with step halving (LDAK sumfuns.c, solve_sums), from theta = start. Optional
// delete-one-group jackknife refits (groups: 0-based, -1 = none).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "la_dispatch.h"
#include "resources.h"

// files: eigen.bin of each block; offset: B + 1 row offsets into the panel; freq: allele frequency; A: panel x Q
// annotations (column 1 = all ones); weight_mode 0 = all 1, 1 = thin, 2 = w_user; reg: regression SNP (0/1).
// [[Rcpp::export]]
Rcpp::List ldak_tagging_cpp(Rcpp::CharacterVector files, Rcpp::IntegerVector offset,
                            const Eigen::Map<Eigen::VectorXd> freq, const Eigen::Map<Eigen::MatrixXd> A, double alpha,
                            int weight_mode, const Eigen::Map<Eigen::VectorXd> w_user, double thin_r2, double n_ref,
                            const Rcpp::IntegerVector& reg, bool normalize, double thresh, int chunk, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int nb = files.size(), M = static_cast<int>(A.rows()), Q = static_cast<int>(A.cols());
  if (offset.size() != nb + 1 || offset[0] != 0 || offset[nb] != M) Rcpp::stop("offset does not match the panel");
  if (freq.size() != M || reg.size() != M) Rcpp::stop("freq / reg do not match the panel");
  const std::vector<std::string> fv = Rcpp::as<std::vector<std::string>>(files);
  const std::vector<int> off(offset.begin(), offset.end());
  const std::vector<int> regv(reg.begin(), reg.end());
  const double c = n_ref > 2 ? 1.0 / (n_ref - 2.0) : 0.0;
  chunk = std::max(1, chunk);
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(M, Q);
  Eigen::VectorXd tg = Eigen::VectorXd::Zero(M), w(M), ex(M);
  std::vector<std::string> err(nb);
  std::vector<int> order(nb);
  for (int i = 0; i < nb; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return off[a + 1] - off[a] > off[b + 1] - off[b]; });
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < nb; ++o) {
    const int i = order[o];
    try {
      const plink4r::EigenBlock ld = plink4r::read_eigen_block(fv[i], thresh);
      const int m = off[i + 1] - off[i], s0 = off[i];
      if (ld.m != m) throw std::runtime_error("the block has " + std::to_string(ld.m) + " eigen rows, the panel " + std::to_string(m));
      const Eigen::MatrixXf& U = ld.U;                                 // m x k
      const int k = static_cast<int>(U.cols());
      const Eigen::VectorXf lam = ld.lambda.cwiseMax(0.0f);
      const Eigen::MatrixXf W = U * lam.asDiagonal();                 // rows w_j = diag(lam) u_j
      const Eigen::MatrixXf Ut = U.transpose();                       // k x m
      const Eigen::VectorXf d = (W.cwiseProduct(U)).rowwise().sum();  // stored diag R_jj
      Eigen::VectorXf sd = Eigen::VectorXf::Ones(m);
      if (normalize) sd = d.cwiseMax(1e-12f).cwiseSqrt().cwiseInverse();
      const ldscr::LaFns& la = ldscr::la();
      auto rows_r2 = [&](int j0, int nr) {
        Eigen::MatrixXf R(nr, m);                                     // R = W[j0:j0+nr, ] U'
        la.sgemm_nn(W.data() + j0, nr, k, m, Ut.data(), m, k, R.data(), nr);
        R = sd.segment(j0, nr).asDiagonal() * R * sd.asDiagonal();
        R = R.array().square().matrix();
        return R;
      };
      // weights
      std::vector<double> wb(m, 1.0);
      if (weight_mode == 2) for (int j = 0; j < m; ++j) wb[j] = w_user[s0 + j];
      if (weight_mode == 1) {
        std::vector<char> keep(m, 1);
        for (int j0 = 0; j0 < m; j0 += chunk) {
          const int nr = std::min(chunk, m - j0);
          const Eigen::MatrixXf r2 = rows_r2(j0, nr);
          for (int a = 0; a < nr; ++a) {
            const int j = j0 + a;
            for (int k = j + 1; k < m; ++k)
              if (keep[k] && r2(a, k) >= thin_r2) { keep[j] = 0; break; }
          }
        }
        for (int j = 0; j < m; ++j) wb[j] = keep[j];
      }
      // expected heritabilities q (m x Q) and the regression indicator
      Eigen::MatrixXf P(m, Q + 1);
      for (int j = 0; j < m; ++j) {
        const double f = freq[s0 + j], het = 2.0 * f * (1.0 - f);
        const double e = het > 0 ? wb[j] * std::pow(het, 1.0 + alpha) : 0.0;
        w[s0 + j] = wb[j];
        ex[s0 + j] = e;
        for (int q = 0; q < Q; ++q) P(j, q) = A(s0 + j, q) * e;
        P(j, Q) = regv[s0 + j];
      }
      const Eigen::RowVectorXd tot = P.cast<double>().colwise().sum();
      for (int j0 = 0; j0 < m; j0 += chunk) {
        const int nr = std::min(chunk, m - j0);
        Eigen::MatrixXf r2 = rows_r2(j0, nr);
        for (int a = 0; a < nr; ++a) r2(a, j0 + a) = 0;               // diagonal handled below
        Eigen::MatrixXf S(nr, Q + 1);                                 // sum_{k != j} r^2 P_k
        la.sgemm_nn(r2.data(), nr, m, nr, P.data(), Q + 1, m, S.data(), nr);
        for (int a = 0; a < nr; ++a) {
          const int j = j0 + a;
          for (int q = 0; q <= Q; ++q) {
            const double pj = P(j, q), val = pj + (1.0 + c) * S(a, q) - c * (tot[q] - pj);
            if (q < Q) V(s0 + j, q) = val; else tg[s0 + j] = val;
          }
        }
      }
    } catch (const std::exception& e) {
      err[i] = e.what();
    }
  }
  for (int i = 0; i < nb; ++i)
    if (!err[i].empty()) Rcpp::stop("%s: %s", fv[i].c_str(), err[i].c_str());
  return Rcpp::List::create(Rcpp::Named("V") = Rcpp::wrap(V), Rcpp::Named("tagging") = Rcpp::wrap(tg),
                            Rcpp::Named("weight") = Rcpp::wrap(w), Rcpp::Named("exp") = Rcpp::wrap(ex));
}

namespace {
struct ShFit {
  Eigen::VectorXd theta;
  double like = 0;
  Eigen::MatrixXd AI;
  bool converged = false;
};

// Per-SNP pieces of the likelihood at expectations E: gradient and observed information.
inline void info_rows(const Eigen::VectorXd& S, const Eigen::MatrixXd& Tm, const Eigen::VectorXd& iw,
                      const Eigen::VectorXd& E, const std::vector<char>& use, Eigen::VectorXd& g, Eigen::MatrixXd& AI) {
  const int n = static_cast<int>(S.size()), K = static_cast<int>(Tm.cols());
  g.setZero(K);
  AI.setZero(K, K);
  for (int j = 0; j < n; ++j) {
    if (!use[j]) continue;
    const double e = E[j], gj = 0.5 * (S[j] - e) * iw[j] / (e * e), hj = (S[j] - 0.5 * e) * iw[j] / (e * e * e);
    const Eigen::VectorXd t = Tm.row(j).transpose();
    g.noalias() += gj * t;
    AI.noalias() += hj * t * t.transpose();
  }
}

inline Eigen::VectorXd expect_of(const Eigen::MatrixXd& Tm, const Eigen::VectorXd& th) {
  Eigen::VectorXd E = (Tm * th).array() + 1.0;
  return E.cwiseMax(1e-6);
}

// Maximises l(theta) - theta' P theta (P = 0: plain SumHer) over the rows with use[j] == 1, by Newton-Raphson
// with step halving on the penalised objective.
ShFit sumher_solve(const Eigen::VectorXd& S, const Eigen::MatrixXd& Tm, const Eigen::VectorXd& iw,
                   const std::vector<char>& use, Eigen::VectorXd theta, const Eigen::MatrixXd& P, double tol,
                   int maxiter) {
  const int n = static_cast<int>(S.size());
  auto objective = [&](const Eigen::VectorXd& E, const Eigen::VectorXd& th) {
    double l = 0;
    for (int j = 0; j < n; ++j) if (use[j]) l += (S[j] / E[j] + std::log(S[j] * E[j])) * iw[j];
    return -0.5 * l - th.dot(P * th);
  };
  ShFit f;
  Eigen::VectorXd E = expect_of(Tm, theta), g;
  Eigen::MatrixXd AI;
  double like = objective(E, theta);
  for (int it = 0; it < maxiter; ++it) {
    info_rows(S, Tm, iw, E, use, g, AI);
    const Eigen::VectorXd step = (AI + 2.0 * P).ldlt().solve(g - 2.0 * P * theta);
    double relax = 1;
    bool moved = false;
    Eigen::VectorXd th2, E2;
    double l2 = like;
    while (relax > 1e-3) {
      th2 = theta + relax * step;
      E2 = expect_of(Tm, th2);
      l2 = objective(E2, th2);
      if (l2 > like - tol) { moved = true; break; }
      relax *= 0.5;
    }
    if (!moved) break;
    const double diff = l2 - like;
    theta = th2; E = E2; like = l2;
    if (std::abs(diff) < tol) { f.converged = true; break; }
  }
  info_rows(S, Tm, iw, E, use, g, AI);
  f.theta = theta; f.like = like; f.AI = AI;
  return f;
}

double median_of(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
}  // namespace

// S: chi-square statistics; Tm: n x K design (n_j / s v_jc, and n_j / s for the intercept); stag: t_j (>= 1);
// group: 0-based SNP groups (contiguous blocks) for the cross-validation and the jackknife;
// pen_d: penalty weight of each coefficient (variance of its annotation column over the panel; 0 = free);
// lambda: NaN = chosen over lambda_grid by leave-one-group-out cross-validation. No ridge when all pen_d are 0.
// Ridge, as in ldsc.strat: penalty lambda * u * sum_c pen_d_c theta_c^2, u = median_c (AI_cc / pen_d_c) at the
// unpenalised fit. Cross-validation criterion: minus the held-out log-likelihood of each group, with the
// leave-one-group-out fit and the held-out likelihood both taken from the second-order expansion at the
// unpenalised fit (one Newton step), so a grid point costs B solves of K x K systems.
// [[Rcpp::export]]
Rcpp::List sumher_cpp(const Eigen::Map<Eigen::VectorXd> S, const Eigen::Map<Eigen::MatrixXd> Tm,
                      const Eigen::Map<Eigen::VectorXd> stag, const Rcpp::IntegerVector& group,
                      const Eigen::Map<Eigen::VectorXd> pen_d, double lambda, const Rcpp::NumericVector& lambda_grid,
                      bool jackknife, double tol, int maxiter) {
  const int n = static_cast<int>(S.size()), K = static_cast<int>(Tm.cols());
  const Eigen::VectorXd S_ = S, iw = stag.cwiseInverse(), pd = pen_d;
  const Eigen::MatrixXd T_ = Tm;
  std::vector<char> all(n, 1);
  const Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(K, K);
  ShFit full = sumher_solve(S_, T_, iw, all, Eigen::VectorXd::Zero(K), Z, tol, maxiter);
  const bool ridge = pd.maxCoeff() > 0;
  const bool have_groups = group.size() == n;
  int B = 0;
  if (have_groups) B = *std::max_element(group.begin(), group.end()) + 1;
  Eigen::MatrixXd P = Z;
  double unit = NA_REAL, lambda_used = NA_REAL;
  Rcpp::NumericVector cv;
  if (ridge) {
    std::vector<double> ratio;
    for (int c = 0; c < K; ++c) if (pd[c] > 0) ratio.push_back(full.AI(c, c) / pd[c]);
    unit = median_of(ratio);
    if (std::isnan(lambda)) {
      if (!have_groups) Rcpp::stop("cross-validation needs SNP groups");
      const Eigen::VectorXd E0 = expect_of(T_, full.theta);
      std::vector<Eigen::VectorXd> gb(B, Eigen::VectorXd::Zero(K));
      std::vector<Eigen::MatrixXd> Ab(B, Eigen::MatrixXd::Zero(K, K));
      for (int j = 0; j < n; ++j) {
        const double e = E0[j], gj = 0.5 * (S_[j] - e) * iw[j] / (e * e), hj = (S_[j] - 0.5 * e) * iw[j] / (e * e * e);
        const Eigen::VectorXd t = T_.row(j).transpose();
        gb[group[j]].noalias() += gj * t;
        Ab[group[j]].noalias() += hj * t * t.transpose();
      }
      Eigen::VectorXd gall = Eigen::VectorXd::Zero(K);
      Eigen::MatrixXd Aall = Eigen::MatrixXd::Zero(K, K);
      for (int b = 0; b < B; ++b) { gall += gb[b]; Aall += Ab[b]; }
      const int G = lambda_grid.size();
      cv = Rcpp::NumericVector(G);
      for (int gi = 0; gi < G; ++gi) {
        const Eigen::MatrixXd Pg = (lambda_grid[gi] * unit * pd).asDiagonal();
        double err = 0;
        for (int b = 0; b < B; ++b) {
          const Eigen::VectorXd d = (Aall - Ab[b] + 2.0 * Pg).ldlt().solve(gall - gb[b] - 2.0 * Pg * full.theta);
          err -= gb[b].dot(d) - 0.5 * d.dot(Ab[b] * d);
        }
        cv[gi] = err;
      }
      const int best = static_cast<int>(std::min_element(cv.begin(), cv.end()) - cv.begin());
      lambda_used = lambda_grid[best];
      cv = cv - cv[best];
    } else {
      lambda_used = lambda;
    }
    P = (lambda_used * unit * pd).asDiagonal();
    if (lambda_used > 0) full = sumher_solve(S_, T_, iw, all, full.theta, P, tol, maxiter);
  }
  Eigen::MatrixXd jk(0, K);
  if (jackknife && have_groups) {
    jk.resize(B, K);
    for (int b = 0; b < B; ++b) {
      std::vector<char> use(n);
      for (int j = 0; j < n; ++j) use[j] = group[j] != b;
      jk.row(b) = sumher_solve(S_, T_, iw, use, full.theta, P, tol, maxiter).theta.transpose();
    }
  }
  // information-based covariance at fixed lambda: H^-1 AI H^-1 with H = AI + 2P (= AI^-1 without ridge)
  const Eigen::MatrixXd Hinv = (full.AI + 2.0 * P).inverse();
  const Eigen::MatrixXd Vt = Hinv * full.AI * Hinv;
  return Rcpp::List::create(Rcpp::Named("theta") = Rcpp::wrap(full.theta), Rcpp::Named("loglik") = full.like,
                            Rcpp::Named("V") = Rcpp::wrap(Vt), Rcpp::Named("converged") = full.converged,
                            Rcpp::Named("jack") = Rcpp::wrap(jk), Rcpp::Named("lambda") = lambda_used,
                            Rcpp::Named("unit") = unit, Rcpp::Named("cv") = cv);
}
