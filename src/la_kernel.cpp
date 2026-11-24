// Matrix-product kernel compiled once per instruction set (src/Makevars) and chosen at run time (la_dispatch.h).
// Adapted from plink4r's la_kernel.cpp: each copy renames the Eigen namespace to Eigen_<isa>, so the copies of
// every Eigen template are distinct symbols and the linker never mixes ISAs. Only Eigen and plain arrays here.
#ifndef LDSCR_LA_ISA
#error "LDSCR_LA_ISA must be defined (generic, avx2 or avx512)"
#endif
#define LDSCR_CAT2(a, b) a##b
#define LDSCR_CAT(a, b) LDSCR_CAT2(a, b)
#define Eigen LDSCR_CAT(Eigen_, LDSCR_LA_ISA)
#include <Eigen/Dense>
#include "la_dispatch.h"

namespace {
using MatF = Eigen::MatrixXf;
void sgemm_nn(const float* A, int M, int K, int lda, const float* B, int N, int ldb, float* Out, int ldo) {
  Eigen::Map<const MatF, 0, Eigen::OuterStride<>> a(A, M, K, Eigen::OuterStride<>(lda));
  Eigen::Map<const MatF, 0, Eigen::OuterStride<>> b(B, K, N, Eigen::OuterStride<>(ldb));
  Eigen::Map<MatF, 0, Eigen::OuterStride<>> o(Out, M, N, Eigen::OuterStride<>(ldo));
  o.noalias() = a * b;
}
}  // namespace

#define LDSCR_LA_TABLE LDSCR_CAT(ldscr_la_, LDSCR_LA_ISA)
#define LDSCR_STR2(x) #x
#define LDSCR_STR(x) LDSCR_STR2(x)
extern const ldscr::LaFns LDSCR_LA_TABLE;
const ldscr::LaFns LDSCR_LA_TABLE = {sgemm_nn, LDSCR_STR(LDSCR_LA_ISA)};
