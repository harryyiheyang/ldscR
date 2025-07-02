#include <RcppArmadillo.h>
#include <random>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List ldsc_bootstrap_cpp(NumericMatrix Zmat, NumericMatrix Nmat, NumericVector LDSC,
                        IntegerVector Block, NumericMatrix GCovInit, NumericMatrix ECovInit,
                        int sampling_time, double z_thresh, double cov_thresh, int M,
                        int n_threads = -1){

  int p = Zmat.ncol();

  // Set number of threads (total_cores/2)
  int num_threads;
  if (n_threads < 1) {
    int max_cores = omp_get_max_threads();
    num_threads = std::max(1, max_cores / 2);  // 默认：最多一半线程
  } else {
    num_threads = n_threads;
  }
  omp_set_num_threads(num_threads);
  Rcout << "Using " << num_threads << " threads for bootstrap\n";

  // Find unique blocks
  arma::ivec blocks = arma::ivec(Block.begin(), Block.size(), false);
  arma::uvec unique_blocks = arma::unique(arma::conv_to<arma::uvec>::from(blocks));
  int nblock = unique_blocks.n_elem;
  int nsub = nblock;  // block-wise bootstrap with replacement: same size

  // Pre-allocate result containers
  std::vector<arma::mat> GCovEstt(sampling_time, arma::mat(p, p, arma::fill::zeros));
  std::vector<arma::mat> ECovEstt(sampling_time, arma::mat(p, p, arma::fill::zeros));

  // Convert R matrices to Armadillo (no copy, just wrapping)
  arma::mat Z = arma::mat(Zmat.begin(), Zmat.nrow(), Zmat.ncol(), false);
  arma::mat N = arma::mat(Nmat.begin(), Nmat.nrow(), Nmat.ncol(), false);
  arma::vec lds = arma::vec(LDSC.begin(), LDSC.size(), false);
  arma::mat GCov_init = arma::mat(GCovInit.begin(), GCovInit.nrow(), GCovInit.ncol(), false);
  arma::mat ECov_init = arma::mat(ECovInit.begin(), ECovInit.nrow(), ECovInit.ncol(), false);

  // Pre-compute block membership for efficiency
  std::vector<std::vector<arma::uword>> block_indices(nblock);
  for (arma::uword j = 0; j < blocks.n_elem; ++j) {
    for (int b = 0; b < nblock; ++b) {
      if (blocks(j) == static_cast<int>(unique_blocks(b))) {
        block_indices[b].push_back(j);
        break;
      }
    }
  }

  // Parallel bootstrap loop
#pragma omp parallel
{
  // Thread-local random number generator with unique seed
  unsigned int seed = 123456 + omp_get_thread_num() * 99991;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> block_dist(0, nblock - 1);

  // Thread-local working variables
  arma::uvec selected_blocks(nblock);
  std::vector<arma::uword> selected_indices;
  selected_indices.reserve(blocks.n_elem * 2); // Pre-reserve memory

  arma::mat X_work, XtX_work;
  arma::vec Xty_work, beta_work;

#pragma omp for schedule(dynamic)
  for (int t = 0; t < sampling_time; ++t) {

    // Clear and reset for this iteration
    selected_indices.clear();

    // Sample blocks with thread-local generator
    for (int i = 0; i < nblock; ++i) {
      selected_blocks(i) = block_dist(gen);
    }

    // Find all SNPs in selected blocks
    for (int i = 0; i < nsub; ++i) {
      int block_idx = selected_blocks(i);
      const std::vector<arma::uword>& indices = block_indices[block_idx];
      selected_indices.insert(selected_indices.end(), indices.begin(), indices.end());
    }

    // Convert to arma::uvec
    arma::uvec row_index(selected_indices.size());
    for (size_t i = 0; i < selected_indices.size(); ++i) {
      row_index(i) = selected_indices[i];
    }

    // Subset data
    arma::mat Zsub = Z.rows(row_index);
    arma::mat Nsub = N.rows(row_index);
    arma::vec lsub = lds.elem(row_index);
    int msub = row_index.n_elem;

    // Diagonal elements (heritability)
    for (int i = 0; i < p; ++i) {
      arma::vec z = Zsub.col(i) % Zsub.col(i); // element-wise square
      z = arma::clamp(z, 0.0, z_thresh);

      arma::vec l = lsub % arma::sqrt(Nsub.col(i) / M) % arma::sqrt(Nsub.col(i) / M);

      // Weight calculation with numerical stability
      arma::vec denom = arma::square(1.0 + l * GCov_init(i, i));
      denom = arma::clamp(denom, 1e-10, arma::datum::inf); // prevent division by zero
      arma::vec w = 1.0 / denom;

      // Design matrix
      X_work.set_size(msub, 2);
      X_work.col(0).ones();
      X_work.col(1) = l;

      // Weighted least squares
      XtX_work = X_work.t() * (X_work.each_col() % w) / msub;
      Xty_work = X_work.t() * (z % w) / msub;

      // Solve with numerical stability check
      if (arma::det(XtX_work) > 1e-12) {
        beta_work = arma::solve(XtX_work, Xty_work);
        GCovEstt[t](i, i) = beta_work(1);
        ECovEstt[t](i, i) = beta_work(0);
      } else {
        // Fallback to regularized solution
        beta_work = arma::solve(XtX_work + 1e-6 * arma::eye(2, 2), Xty_work);
        GCovEstt[t](i, i) = beta_work(1);
        ECovEstt[t](i, i) = beta_work(0);
      }
    }

    // Off-diagonal elements (genetic covariances)
    for (int i = 0; i < p - 1; ++i) {
      for (int j = i + 1; j < p; ++j) {
        arma::vec z = Zsub.col(i) % Zsub.col(j); // element-wise product
        z = arma::clamp(z, -cov_thresh, cov_thresh);

        arma::vec l = lsub % arma::sqrt(Nsub.col(i) / M) % arma::sqrt(Nsub.col(j) / M);
        arma::vec li = lsub % (Nsub.col(i) / M);
        arma::vec lj = lsub % (Nsub.col(j) / M);

        // Weight calculation with numerical stability
        arma::vec term1 = (1.0 + li * GCov_init(i, i)) % (1.0 + lj * GCov_init(j, j));
        arma::vec term2 = arma::square(l * GCov_init(i, j) + ECov_init(i, j));
        arma::vec denom = term1 + term2;
        denom = arma::clamp(denom, 1e-10, arma::datum::inf); // prevent division by zero
        arma::vec w = 1.0 / denom;

        // Design matrix
        X_work.set_size(msub, 2);
        X_work.col(0).ones();
        X_work.col(1) = l;

        // Weighted least squares
        XtX_work = X_work.t() * (X_work.each_col() % w) / msub;
        Xty_work = X_work.t() * (z % w) / msub;

        // Solve with numerical stability check
        if (arma::det(XtX_work) > 1e-12) {
          beta_work = arma::solve(XtX_work, Xty_work);
          GCovEstt[t](i, j) = GCovEstt[t](j, i) = beta_work(1);
          ECovEstt[t](i, j) = ECovEstt[t](j, i) = beta_work(0);
        } else {
          // Fallback to regularized solution
          beta_work = arma::solve(XtX_work + 1e-6 * arma::eye(2, 2), Xty_work);
          GCovEstt[t](i, j) = GCovEstt[t](j, i) = beta_work(1);
          ECovEstt[t](i, j) = ECovEstt[t](j, i) = beta_work(0);
        }
      }
    }

    // Progress indicator (thread-safe)
    if ((t + 1) % 50 == 0) {
#pragma omp critical
{
  Rcout << "Completed " << (t + 1) << " / " << sampling_time << " bootstrap samples\n";
}
    }
  } // end of parallel for loop

  // Check for user interrupt in each thread
#pragma omp master
{
  Rcpp::checkUserInterrupt();
}

} // end of parallel region

return List::create(Named("GCovEstt") = GCovEstt,
                    Named("ECovEstt") = ECovEstt);
}
