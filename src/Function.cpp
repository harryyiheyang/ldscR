#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 辅助函数：计算向量与矩阵每一列的点积
// [[Rcpp::export]]
arma::mat vectorMatrixDotProduct(const arma::vec& v, const arma::mat& M) {
  if (v.size() != M.n_rows) {
    stop("The vector length and the number of rows in the matrix must be the same.");
  }

  arma::mat results(M.n_rows, M.n_cols);
  for (unsigned int i = 0; i < M.n_cols; ++i) {
    results.col(i) = v % M.col(i);
  }

  return results;
}


// 辅助函数：计算自乘积和成对乘积
// [[Rcpp::export]]
List selfPairwiseProduct(const mat& M) {
int M_rows = M.n_rows, M_cols = M.n_cols;
mat selfProducts = M % M;  // 每个元素自身的乘积
int numCols = (M_cols * (M_cols - 1)) / 2;
mat pairwiseProducts(M_rows, numCols);
int k = 0;
for (int j = 0; j < M_cols; ++j) {
for (int s = j + 1; s < M_cols; ++s) {
pairwiseProducts.col(k++) = M.col(j) % M.col(s);
}
}
return List::create(Named("selfProducts") = selfProducts,
              Named("pairwiseProducts") = pairwiseProducts);
}


// 标准LDSC函数(无权重)
// [[Rcpp::export]]
List ldsc(const mat& selfProductZ, const mat& pairwiseProductZ,
          const mat& selfProductN, const mat& pairwiseProductN) {
  int p = selfProductZ.n_cols;
  vec Hest1(p), Hest2(p);
  vec CoVest1(pairwiseProductZ.n_cols), CoVest2(pairwiseProductZ.n_cols);

  // 执行回归分析
  for (int i = 0; i < p; ++i) {
    mat X = join_horiz(ones<vec>(selfProductN.n_rows), selfProductN.col(i));
    vec y = selfProductZ.col(i);
    vec beta = solve(X, y);
    Hest1(i) = beta(0);
    Hest2(i) = beta(1);
  }

  for (unsigned int i = 0; i < pairwiseProductZ.n_cols; ++i) {
    mat X = join_horiz(ones<vec>(pairwiseProductN.n_rows), pairwiseProductN.col(i));
    vec y = pairwiseProductZ.col(i);
    vec beta = solve(X, y);
    CoVest1(i) = beta(0);
    CoVest2(i) = beta(1);
  }

  // 构建GCovEst1和ECovEst1矩阵
  mat GCovEst1 = eye<mat>(p, p);
  mat ECovEst1 = eye<mat>(p, p);
  GCovEst1.diag() = Hest2;
  ECovEst1.diag() = Hest1;

  if (p > 1) {
    vec CoVest1_reversed = reverse(CoVest1);
    vec CoVest2_reversed = reverse(CoVest2);
    int k = 0;
    for (int i = 1; i < p; ++i) {
      for (int j = 0; j < i; ++j) {
        GCovEst1(i, j) = CoVest2_reversed(k);
        ECovEst1(i, j) = CoVest1_reversed(k);
        ++k;
      }
    }
  }

  return List::create(Named("GCovEst1") = GCovEst1, Named("ECovEst1") = ECovEst1);
}

// 加权LDSC函数
// [[Rcpp::export]]
List weightedLdsc(const mat& selfProductZ, const mat& pairwiseProductZ,
                  const mat& selfProductN, const mat& pairwiseProductN,
                  const mat& GCovEst, const mat& ECovEst) {
  int p = selfProductZ.n_cols;
  vec Hest1(p), Hest2(p);
  vec CoVest1(pairwiseProductZ.n_cols), CoVest2(pairwiseProductZ.n_cols);

  // 执行加权回归分析
  for (int i = 0; i < p; ++i) {
    mat X = join_horiz(ones<vec>(selfProductN.n_rows), selfProductN.col(i));
    vec y = selfProductZ.col(i);
    vec w = square(1 + selfProductN.col(i) * GCovEst(i,i));
    vec w_sqrt = sqrt(1 / w);
    X.each_col() %= w_sqrt;
    y %= w_sqrt;
    vec beta = solve(X, y);
    Hest1(i) = beta(0);
    Hest2(i) = beta(1);
  }

  for (unsigned int i = 0; i < pairwiseProductZ.n_cols; ++i) {
    mat X = join_horiz(ones<vec>(pairwiseProductN.n_rows), pairwiseProductN.col(i));
    vec y = pairwiseProductZ.col(i);
    int j = i + 1;
    vec w = (1 + selfProductN.col(i) * GCovEst(i,i)) % (1 + selfProductN.col(j) * GCovEst(j,j)) + square(pairwiseProductN.col(i) * GCovEst(i,j) + ECovEst(i,j));
    vec w_sqrt = sqrt(1 / w);
    X.each_col() %= w_sqrt;
    y %= w_sqrt;
    vec beta = solve(X, y);
    CoVest1(i) = beta(0);
    CoVest2(i) = beta(1);
  }

  // 构建GCovEst1和ECovEst1矩阵
  mat GCovEst1 = eye<mat>(p, p);
  mat ECovEst1 = eye<mat>(p, p);
  GCovEst1.diag() = Hest2;
  ECovEst1.diag() = Hest1;

  if (p > 1) {
    vec CoVest1_reversed = reverse(CoVest1);
    vec CoVest2_reversed = reverse(CoVest2);
    int k = 0;
    for (int i = 1; i < p; ++i) {
      for (int j = 0; j < i; ++j) {
        GCovEst1(i, j) = CoVest2_reversed(k);
        ECovEst1(i, j) = CoVest1_reversed(k);
        ++k;
      }
    }
  }

  return List::create(Named("GCovEst1") = GCovEst1, Named("ECovEst1") = ECovEst1);
}
