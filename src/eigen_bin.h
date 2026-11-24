// eigen.bin reader (the file of a block is named by R: .ld_block_file). Two layouts, told apart by file size:
//  SBayesRC : int32 m; int32 k; float32 sumLambda; float32 thresh; float32 lam[k]; float32 U[m*k] (col-major)
//  plink4r  : int32 m; int32 k; float64 sumLambda; float32 thresh; float64 lam[k]; float32 scale[k]; int16 U[m*k] (col-major, u = q*scale)
#ifndef PLINK4R_EIGEN_BIN_H
#define PLINK4R_EIGEN_BIN_H

#include <RcppEigen.h>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace plink4r {

struct EigenBlock {
  int m = 0, k = 0;
  double sum_lambda = 0;
  float file_thresh = 0;
  Eigen::VectorXf lambda;
  Eigen::MatrixXf U;
};

template <typename T>
inline void read_n(std::ifstream& in, T* dst, std::size_t n, const std::string& path) {
  in.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(sizeof(T) * n));
  if (!in) throw std::runtime_error("short read: " + path);
}

// thresh <= 0 : keep every stored eigenpair; otherwise keep the smallest k reaching thresh * sumLambda.
inline EigenBlock read_eigen_block(const std::string& path, double thresh) {
  std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("cannot open " + path);
  const std::uint64_t size = static_cast<std::uint64_t>(in.tellg());
  in.seekg(0);
  int32_t m = 0, kf = 0;
  read_n(in, &m, 1, path);
  read_n(in, &kf, 1, path);
  if (m <= 0 || kf <= 0 || kf > m) throw std::runtime_error("bad dimensions: " + path);
  const std::uint64_t sz_sbrc = 16ull + 4ull * kf + 4ull * m * kf;
  const std::uint64_t sz_p4r = 20ull + 12ull * kf + 2ull * m * kf;
  EigenBlock b;
  b.m = m;
  std::vector<double> lam(kf);
  std::vector<float> scale;
  bool p4r;
  if (size == sz_p4r) {
    p4r = true;
    read_n(in, &b.sum_lambda, 1, path);
    read_n(in, &b.file_thresh, 1, path);
    read_n(in, lam.data(), kf, path);
    scale.resize(kf);
    read_n(in, scale.data(), kf, path);
  } else if (size == sz_sbrc) {
    p4r = false;
    float s = 0;
    read_n(in, &s, 1, path);
    b.sum_lambda = s;
    read_n(in, &b.file_thresh, 1, path);
    std::vector<float> lf(kf);
    read_n(in, lf.data(), kf, path);
    for (int j = 0; j < kf; ++j) lam[j] = lf[j];
  } else {
    throw std::runtime_error("unknown eigen.bin layout (size mismatch): " + path);
  }
  int k = kf;
  if (thresh > 0) {
    if (thresh > b.file_thresh + 1e-6) throw std::runtime_error("thresh exceeds the file cutoff: " + path);
    double cum = 0;
    for (int j = 0; j < kf; ++j) {
      cum += lam[j];
      if (cum >= thresh * b.sum_lambda) { k = j + 1; break; }
    }
  }
  b.k = k;
  b.lambda.resize(k);
  for (int j = 0; j < k; ++j) b.lambda[j] = static_cast<float>(lam[j]);
  b.U.resize(m, k);
  if (p4r) {
    std::vector<int16_t> q(static_cast<std::size_t>(m) * k);
    read_n(in, q.data(), q.size(), path);
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < m; ++i) b.U(i, j) = q[static_cast<std::size_t>(j) * m + i] * scale[j];
  } else {
    read_n(in, b.U.data(), static_cast<std::size_t>(m) * k, path);
  }
  return b;
}

}  // namespace plink4r
#endif
